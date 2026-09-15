import assert from "node:assert/strict";
import test from "node:test";

import { buildStreamlines, createWindSampler, wrapLongitude } from "./static/streamlines.js";

const grid = {
  mu: [-0.75, 0, 0.75],
  nlon: [4, 4, 4],
  ring_offsets: [0, 4, 8, 12],
  point_count: 12,
};

test("wrapLongitude keeps angles in the map interval", () => {
  assert.ok(Math.abs(wrapLongitude(3 * Math.PI) + Math.PI) < 1e-12);
  assert.ok(Math.abs(wrapLongitude(-2.5 * Math.PI) + Math.PI / 2) < 1e-12);
});

test("wind sampler interpolates between rings and across the date line", () => {
  const u = Float32Array.from([
    1, 1, 1, 1,
    3, 3, 3, 3,
    5, 5, 5, 5,
  ]);
  const v = new Float32Array(12);
  const sample = createWindSampler(grid, u, v);
  assert.equal(sample(Math.PI - 1e-8, 0).u, 3);
  assert.equal(sample(-Math.PI + 1e-8, 0).u, 3);
  assert.ok(Math.abs(sample(0, Math.asin(0.375)).u - 4) < 1e-6);
});

test("uniform eastward flow produces zonal streamlines", () => {
  const u = new Float32Array(12).fill(20);
  const v = new Float32Array(12);
  const result = buildStreamlines(grid, u, v, "low");
  assert.equal(result.maximumSpeed, 20);
  assert.ok(result.lines.length > 0);
  const middle = result.lines[Math.floor(result.lines.length / 2)];
  assert.ok(middle.length >= 8);
  const latitude = middle[0].latitude;
  for (const point of middle) assert.ok(Math.abs(point.latitude - latitude) < 1e-10);
  assert.ok(middle.at(-1).longitude > middle[0].longitude);
});
