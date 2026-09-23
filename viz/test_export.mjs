import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import test from "node:test";
import { crc32, createZip } from "./static/export.js";

test("crc32 matches the reference value", () => {
  assert.equal(crc32(new TextEncoder().encode("123456789")), 0xcbf43926);
});

test("createZip writes an archive that Python's zipfile reads back", async () => {
  const csv = `a,b\n${Array.from({ length: 200 }, (_, i) => `${i},${i * i}`).join("\n")}\n`;
  const png = new Uint8Array([0x89, 0x50, 0x4e, 0x47, 1, 2, 3, 4]);
  const zip = await createZip([
    { name: "run/README.txt", data: "日本語のテキスト\n" },
    { name: "run/csv/data.csv", data: csv },
    { name: "run/maps/a.png", data: new Blob([png]) },
  ]);
  const path = join(mkdtempSync(join(tmpdir(), "bespoke-zip-")), "out.zip");
  writeFileSync(path, new Uint8Array(await zip.arrayBuffer()));
  const script = `
import json, sys, zipfile
with zipfile.ZipFile(sys.argv[1]) as z:
    assert z.testzip() is None
    print(json.dumps({i.filename: [i.compress_type, z.read(i).decode("latin-1")] for i in z.infolist()}))
`;
  const result = JSON.parse(execFileSync("python3", ["-c", script, path], { encoding: "utf8" }));
  assert.deepEqual(Object.keys(result), ["run/README.txt", "run/csv/data.csv", "run/maps/a.png"]);
  assert.equal(Buffer.from(result["run/README.txt"][1], "latin1").toString("utf8"), "日本語のテキスト\n");
  assert.equal(result["run/csv/data.csv"][0], 8);
  assert.equal(result["run/csv/data.csv"][1], csv);
  assert.equal(result["run/maps/a.png"][0], 0);
  assert.deepEqual([...Buffer.from(result["run/maps/a.png"][1], "latin1")], [...png]);
});
