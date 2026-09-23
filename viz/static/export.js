// Export helpers of the land-sea visualizer: a small ZIP writer, and figures
// (equirectangular maps and analysis charts) drawn on a 2D canvas with their
// title and legend so that the saved image can be read on its own.

// ---------------------------------------------------------------------------
// ZIP
// ---------------------------------------------------------------------------

const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n += 1) {
    let c = n;
    for (let k = 0; k < 8; k += 1) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

export function crc32(bytes) {
  let crc = 0xffffffff;
  for (let i = 0; i < bytes.length; i += 1) crc = CRC_TABLE[(crc ^ bytes[i]) & 0xff] ^ (crc >>> 8);
  return (crc ^ 0xffffffff) >>> 0;
}

async function deflateRaw(bytes) {
  const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream("deflate-raw"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function toBytes(data) {
  if (data instanceof Uint8Array) return data;
  if (typeof data === "string") return new TextEncoder().encode(data);
  if (data instanceof Blob) return new Uint8Array(await data.arrayBuffer());
  if (data instanceof ArrayBuffer) return new Uint8Array(data);
  throw new TypeError("ZIP に入れるデータは文字列・Blob・バイト列のいずれかです");
}

/**
 * A ZIP archive (Blob) of files [{name, data}], data a string (UTF-8), Blob or
 * bytes. Text is deflated when CompressionStream is available; images, already
 * compressed, are stored.
 */
export async function createZip(files, date = new Date()) {
  const encoder = new TextEncoder();
  const time = (date.getHours() << 11) | (date.getMinutes() << 5) | Math.floor(date.getSeconds() / 2);
  const day = ((Math.max(1980, date.getFullYear()) - 1980) << 9) | ((date.getMonth() + 1) << 5) | date.getDate();
  const canDeflate = typeof CompressionStream === "function";
  const parts = [];
  const central = [];
  let offset = 0;
  for (const file of files) {
    const name = encoder.encode(file.name);
    const raw = await toBytes(file.data);
    const crc = crc32(raw);
    let method = 0;
    let body = raw;
    if (canDeflate && !/\.(png|jpe?g|zip|gz)$/i.test(file.name) && raw.length > 64) {
      const packed = await deflateRaw(raw);
      if (packed.length < raw.length) { method = 8; body = packed; }
    }
    const local = new DataView(new ArrayBuffer(30));
    local.setUint32(0, 0x04034b50, true);
    local.setUint16(4, 20, true);
    local.setUint16(6, 0x0800, true); // UTF-8 names
    local.setUint16(8, method, true);
    local.setUint16(10, time, true);
    local.setUint16(12, day, true);
    local.setUint32(14, crc, true);
    local.setUint32(18, body.length, true);
    local.setUint32(22, raw.length, true);
    local.setUint16(26, name.length, true);
    local.setUint16(28, 0, true);
    parts.push(local, name, body);
    const entry = new DataView(new ArrayBuffer(46));
    entry.setUint32(0, 0x02014b50, true);
    entry.setUint16(4, 20, true);
    entry.setUint16(6, 20, true);
    entry.setUint16(8, 0x0800, true);
    entry.setUint16(10, method, true);
    entry.setUint16(12, time, true);
    entry.setUint16(14, day, true);
    entry.setUint32(16, crc, true);
    entry.setUint32(20, body.length, true);
    entry.setUint32(24, raw.length, true);
    entry.setUint16(28, name.length, true);
    entry.setUint32(42, offset, true);
    central.push(entry, name);
    offset += 30 + name.length + body.length;
  }
  const centralSize = central.reduce((sum, part) => sum + part.byteLength, 0);
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054b50, true);
  end.setUint16(8, files.length, true);
  end.setUint16(10, files.length, true);
  end.setUint32(12, centralSize, true);
  end.setUint32(16, offset, true);
  return new Blob([...parts, ...central, end], { type: "application/zip" });
}

// ---------------------------------------------------------------------------
// Figures
// ---------------------------------------------------------------------------

export const FIGURE = {
  background: "#0c1b27",
  text: "#dfe7ec",
  muted: "#8b9aa6",
  line: "rgba(174, 198, 214, 0.32)",
  font: "\"Hiragino Sans\", \"Noto Sans JP\", system-ui, sans-serif",
  mono: "ui-monospace, SFMono-Regular, Menlo, \"Hiragino Sans\", \"Noto Sans JP\", monospace",
  scale: 2,
};

const font = (size, weight = 400) => `${weight} ${size}px ${FIGURE.font}`;

/** Lines of text that fit in width (breaks at spaces and between CJK characters). */
function wrapText(ctx, text, width) {
  if (!text) return [];
  const tokens = text.match(/[　-鿿＀-￯][、。・）」]?|[^\s　-鿿＀-￯]+\s*|\s+/g) ?? [text];
  const lines = [];
  let line = "";
  for (const token of tokens) {
    if (line && ctx.measureText(line + token).width > width) {
      lines.push(line.trimEnd());
      line = token.trimStart();
    } else {
      line += token;
    }
  }
  if (line.trim()) lines.push(line.trimEnd());
  return lines;
}

/**
 * Legend layout. legend: { kind: "categorical", entries: [{color, label}] } or
 * { kind: "continuous", stops, min, max, unit, ticks, missing: label|null }.
 * Returns { height, draw(ctx, x, y) } for the given width.
 */
function layoutLegend(ctx, legend, width) {
  if (!legend) return { height: 0, draw() {} };
  if (legend.kind === "categorical") {
    ctx.font = font(12);
    const rowHeight = 20;
    const items = [];
    let x = 0;
    let row = 0;
    for (const entry of legend.entries) {
      const itemWidth = 18 + ctx.measureText(entry.label).width;
      if (x > 0 && x + itemWidth > width) { x = 0; row += 1; }
      items.push({ ...entry, x, row });
      x += itemWidth + 18;
    }
    return {
      height: (row + 1) * rowHeight,
      draw(target, left, top) {
        target.font = font(12);
        target.textBaseline = "middle";
        target.textAlign = "left";
        for (const item of items) {
          const y = top + item.row * rowHeight + rowHeight / 2;
          if (item.line) {
            // Series key of a line chart (dashed for reference series).
            target.strokeStyle = item.color;
            target.lineWidth = 2.5;
            target.setLineDash(item.dash ? [3, 2] : []);
            target.beginPath();
            target.moveTo(left + item.x, y);
            target.lineTo(left + item.x + 12, y);
            target.stroke();
            target.setLineDash([]);
          } else {
            target.fillStyle = item.color;
            target.fillRect(left + item.x, y - 6, 12, 12);
            target.strokeStyle = "rgba(0, 0, 0, 0.35)";
            target.lineWidth = 1;
            target.strokeRect(left + item.x + 0.5, y - 5.5, 11, 11);
          }
          target.fillStyle = FIGURE.text;
          target.fillText(item.label, left + item.x + 18, y);
        }
      },
    };
  }
  const barWidth = Math.min(520, width - 160);
  return {
    height: 44,
    draw(target, left, top) {
      const gradient = target.createLinearGradient(left, 0, left + barWidth, 0);
      legend.stops.forEach((color, i) => gradient.addColorStop(i / (legend.stops.length - 1), color));
      target.fillStyle = gradient;
      target.fillRect(left, top + 4, barWidth, 14);
      target.strokeStyle = FIGURE.line;
      target.strokeRect(left + 0.5, top + 4.5, barWidth - 1, 13);
      target.font = font(11);
      target.fillStyle = FIGURE.text;
      target.textBaseline = "top";
      const span = legend.max - legend.min;
      for (const tick of legend.ticks) {
        const x = left + ((tick.value - legend.min) / span) * barWidth;
        target.fillStyle = FIGURE.muted;
        target.fillRect(Math.round(x), top + 18, 1, 4);
        target.fillStyle = FIGURE.text;
        target.textAlign = "center";
        target.fillText(tick.label, x, top + 24);
      }
      target.textAlign = "left";
      target.textBaseline = "middle";
      let x = left + barWidth + 12;
      if (legend.unit && legend.unit !== "1") {
        target.fillText(legend.unit, x, top + 11);
        x += target.measureText(legend.unit).width + 24;
      }
      if (legend.missing) {
        target.fillStyle = legend.missingColor;
        target.fillRect(x, top + 5, 12, 12);
        target.fillStyle = FIGURE.text;
        target.fillText(legend.missing, x + 18, top + 11);
      }
    },
  };
}

/**
 * A figure: title and subtitle, a body of bodyWidth × bodyHeight drawn by
 * drawBody(ctx, x, y), and the legend below. Returns a canvas.
 */
export function composeFigure({ title, subtitle = "", legend = null, bodyWidth, bodyHeight, drawBody, footer = "" }) {
  const pad = 20;
  const width = bodyWidth + 2 * pad;
  const measure = document.createElement("canvas").getContext("2d");
  measure.font = font(16, 600);
  const titleLines = wrapText(measure, title, width - 2 * pad);
  measure.font = font(12);
  const subtitleLines = wrapText(measure, subtitle, width - 2 * pad);
  const footerLines = wrapText(measure, footer, width - 2 * pad);
  const legendLayout = layoutLegend(measure, legend, width - 2 * pad);
  const headerHeight = titleLines.length * 22 + subtitleLines.length * 17 + 10;
  const legendHeight = legendLayout.height ? legendLayout.height + 14 : 0;
  const footerHeight = footerLines.length ? footerLines.length * 16 + 6 : 0;
  const height = pad + headerHeight + bodyHeight + legendHeight + footerHeight + pad;
  const canvas = document.createElement("canvas");
  canvas.width = Math.round(width * FIGURE.scale);
  canvas.height = Math.round(height * FIGURE.scale);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(FIGURE.scale, 0, 0, FIGURE.scale, 0, 0);
  ctx.fillStyle = FIGURE.background;
  ctx.fillRect(0, 0, width, height);
  let y = pad;
  ctx.textAlign = "left";
  ctx.textBaseline = "top";
  ctx.fillStyle = FIGURE.text;
  ctx.font = font(16, 600);
  for (const line of titleLines) { ctx.fillText(line, pad, y); y += 22; }
  ctx.fillStyle = FIGURE.muted;
  ctx.font = font(12);
  for (const line of subtitleLines) { ctx.fillText(line, pad, y); y += 17; }
  y += 10;
  ctx.save();
  drawBody(ctx, pad, y);
  ctx.restore();
  y += bodyHeight;
  if (legendHeight) {
    y += 14;
    legendLayout.draw(ctx, pad, y);
    y += legendLayout.height;
  }
  if (footerLines.length) {
    y += 6;
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.fillStyle = FIGURE.muted;
    ctx.font = font(11);
    for (const line of footerLines) { ctx.fillText(line, pad, y); y += 16; }
  }
  return canvas;
}

const LON_TICKS = [
  [-180, "180°"], [-120, "西経120°"], [-60, "西経60°"], [0, "0°"], [60, "東経60°"], [120, "東経120°"], [180, "180°"],
];
const LAT_TICKS = [[60, "北緯60°"], [30, "北緯30°"], [0, "赤道"], [-30, "南緯30°"], [-60, "南緯60°"]];

/**
 * Equirectangular map of the native cells (prime meridian at the centre).
 * colors: Float32Array of sRGB triplets in [0, 1] per cell. coastline: segments
 * [lon0, lat0, lon1, lat1] in map longitude; sites: [{index, label, custom}].
 */
export function mapFigure({ grid, colors, coastline = [], sites = [], title, subtitle, legend, footer, width = 960 }) {
  const axisLeft = 58;
  const axisBottom = 22;
  const mapWidth = width - axisLeft - 8;
  const mapHeight = mapWidth / 2;
  return composeFigure({
    title, subtitle, legend, footer, bodyWidth: width, bodyHeight: mapHeight + axisBottom,
    drawBody(ctx, left, top) {
      const x0 = left + axisLeft;
      const px = (lon) => x0 + ((lon + 180) / 360) * mapWidth;
      const py = (lat) => top + ((90 - lat) / 180) * mapHeight;
      // Cells are painted on an unscaled canvas so that neighbours share pixel edges.
      const s = FIGURE.scale;
      const cells = document.createElement("canvas");
      cells.width = Math.round(mapWidth * s);
      cells.height = Math.round(mapHeight * s);
      const cctx = cells.getContext("2d");
      const cx = (lon) => Math.round(((lon + 180) / 360) * cells.width);
      const cy = (lat) => Math.round(((90 - lat) / 180) * cells.height);
      const rgb = (i) => `rgb(${Math.round(colors[i * 3] * 255)}, ${Math.round(colors[i * 3 + 1] * 255)}, ${Math.round(colors[i * 3 + 2] * 255)})`;
      for (let j = 0; j < grid.nlat; j += 1) {
        const cellWidth = 360 / grid.nlon[j];
        const yTop = cy(grid.latEdges[j + 1]);
        const yBottom = cy(grid.latEdges[j]);
        for (let k = 0; k < grid.nlon[j]; k += 1) {
          const cell = grid.offsets[j] + k;
          cctx.fillStyle = rgb(cell);
          const west = grid.mapLon[cell] - cellWidth / 2;
          const east = grid.mapLon[cell] + cellWidth / 2;
          const spans = west < -180 ? [[-180, east], [west + 360, 180]]
            : east > 180 ? [[west, 180], [-180, east - 360]] : [[west, east]];
          // One extra pixel to the east hides rounding gaps between neighbours; the next cell paints over it.
          for (const [a, b] of spans) cctx.fillRect(cx(a), yTop, cx(b) - cx(a) + 1, yBottom - yTop);
        }
      }
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(cells, x0, top, mapWidth, mapHeight);
      ctx.save();
      ctx.beginPath();
      ctx.rect(x0, top, mapWidth, mapHeight);
      ctx.clip();
      ctx.strokeStyle = "rgba(255, 255, 255, 0.22)";
      ctx.lineWidth = 0.6;
      ctx.beginPath();
      for (let lon = -120; lon <= 120; lon += 60) { ctx.moveTo(px(lon), top); ctx.lineTo(px(lon), top + mapHeight); }
      for (let lat = -60; lat <= 60; lat += 30) { ctx.moveTo(x0, py(lat)); ctx.lineTo(x0 + mapWidth, py(lat)); }
      ctx.stroke();
      if (coastline.length) {
        ctx.strokeStyle = "rgba(20, 24, 28, 0.95)";
        ctx.lineWidth = 1;
        ctx.lineCap = "round";
        ctx.beginPath();
        for (const [lon0, lat0, lon1, lat1] of coastline) { ctx.moveTo(px(lon0), py(lat0)); ctx.lineTo(px(lon1), py(lat1)); }
        ctx.stroke();
      }
      ctx.restore();
      ctx.strokeStyle = FIGURE.line;
      ctx.lineWidth = 1;
      ctx.strokeRect(x0 + 0.5, top + 0.5, mapWidth - 1, mapHeight - 1);
      for (const site of sites) {
        const x = px(grid.mapLon[site.index]);
        const y = py(grid.lat[site.index]);
        ctx.font = `700 10px ${FIGURE.mono}`;
        const w = Math.max(17, ctx.measureText(site.label).width + 8);
        ctx.fillStyle = site.custom ? "#c8f06b" : "#ffffff";
        ctx.strokeStyle = "#111111";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.roundRect(x - w / 2, y - 8.5, w, 17, 8.5);
        ctx.fill();
        ctx.stroke();
        ctx.fillStyle = "#111111";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(site.label, x, y + 0.5);
      }
      ctx.font = font(11);
      ctx.fillStyle = FIGURE.muted;
      ctx.textBaseline = "top";
      ctx.textAlign = "center";
      for (const [lon, label] of LON_TICKS) ctx.fillText(label, px(lon), top + mapHeight + 6);
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      for (const [lat, label] of LAT_TICKS) ctx.fillText(label, x0 - 6, py(lat));
    },
  });
}

/** A figure around an existing canvas (chart or WebGL view) of css size width × height. */
export function canvasFigure({ source, width, height, title, subtitle, legend, footer }) {
  return composeFigure({
    title, subtitle, legend, footer, bodyWidth: width, bodyHeight: height,
    drawBody(ctx, left, top) { ctx.drawImage(source, left, top, width, height); },
  });
}

export function canvasToBlob(canvas) {
  return new Promise((resolve, reject) => canvas.toBlob((blob) => (blob ? resolve(blob) : reject(new Error("PNG を作れませんでした"))), "image/png"));
}
