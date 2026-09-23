// Small canvas chart kit for the analysis panels: line charts with a hover
// crosshair, climographs, and a filled-contour latitude-pressure section.

export const THEME = {
  surface: "#0c1b27",
  text: "#dfe7ec",
  muted: "#8b9aa6",
  grid: "rgba(174, 198, 214, 0.13)",
  axis: "rgba(174, 198, 214, 0.32)",
  font: "10px ui-monospace, SFMono-Regular, Menlo, \"Hiragino Sans\", \"Noto Sans JP\", monospace",
};
// Categorical order for the dark surface (dataviz reference palette).
export const SERIES = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"];

function niceStep(min, max, count) {
  const raw = (max - min) / Math.max(1, count);
  const power = 10 ** Math.floor(Math.log10(raw));
  return [1, 2, 2.5, 5, 10].map((m) => m * power).find((candidate) => candidate >= raw) ?? 10 * power;
}

/** A range of round numbers that encloses [min, max] (ticks fall on both ends). */
export function niceRange(min, max, count = 5) {
  if (!Number.isFinite(min) || !Number.isFinite(max)) return [0, 1];
  if (min === max) { min -= 1; max += 1; }
  const step = niceStep(min, max, count);
  return [Math.floor(min / step + 1e-9) * step, Math.ceil(max / step - 1e-9) * step];
}

export function niceTicks(min, max, count = 5) {
  if (!Number.isFinite(min) || !Number.isFinite(max)) return [];
  if (min === max) { min -= 1; max += 1; }
  const step = niceStep(min, max, count);
  const ticks = [];
  for (let value = Math.ceil(min / step) * step; value <= max + step * 1e-9; value += step) {
    ticks.push(Math.abs(value) < step * 1e-9 ? 0 : value);
  }
  return ticks;
}

export function formatNumber(value, digits = 4) {
  if (!Number.isFinite(value)) return "—";
  const magnitude = Math.abs(value);
  if (magnitude !== 0 && (magnitude < 1e-3 || magnitude >= 1e6)) return value.toExponential(2);
  return Number(value.toPrecision(digits)).toLocaleString("en-US", { maximumFractionDigits: 6 });
}

export const LATITUDE_TICKS = [-90, -60, -30, 0, 30, 60, 90].map((value) => ({
  value, label: value === 0 ? "赤道" : `${value < 0 ? "南緯" : "北緯"}${Math.abs(value)}°`,
}));

function setupCanvas(canvas, height) {
  const width = Math.max(10, canvas.parentElement.clientWidth);
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.style.height = `${height}px`;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = THEME.surface;
  ctx.fillRect(0, 0, width, height);
  ctx.font = THEME.font;
  return { ctx, width, height };
}

function extent(values, include = []) {
  let min = Infinity;
  let max = -Infinity;
  for (const value of [...values, ...include]) {
    if (!Number.isFinite(value)) continue;
    min = Math.min(min, value);
    max = Math.max(max, value);
  }
  return min <= max ? [min, max] : [0, 1];
}

class BaseChart {
  constructor(container, spec) {
    this.container = container;
    this.spec = spec;
    this.canvas = document.createElement("canvas");
    this.canvas.className = "chart-canvas";
    this.tooltip = document.createElement("div");
    this.tooltip.className = "chart-tooltip";
    this.tooltip.hidden = true;
    container.append(this.canvas, this.tooltip);
    this.hover = null;
    this.canvas.addEventListener("pointermove", (event) => this.onPointer(event));
    this.canvas.addEventListener("pointerleave", () => { this.hover = null; this.tooltip.hidden = true; this.draw(); });
    this.observer = new ResizeObserver(() => this.draw());
    this.observer.observe(container);
  }

  destroy() { this.observer.disconnect(); }

  onPointer(event) {
    const rect = this.canvas.getBoundingClientRect();
    this.hover = { x: event.clientX - rect.left, y: event.clientY - rect.top };
    this.draw();
  }

  showTooltip(lines, x, y) {
    this.tooltip.replaceChildren(...lines.map((line) => {
      const row = document.createElement("div");
      if (line.color) {
        const key = document.createElement("i");
        key.style.background = line.color;
        row.append(key);
      }
      row.append(document.createTextNode(line.text));
      return row;
    }));
    this.tooltip.hidden = false;
    const box = this.container.getBoundingClientRect();
    const left = x + 14 + this.tooltip.offsetWidth > box.width ? x - 14 - this.tooltip.offsetWidth : x + 14;
    this.tooltip.style.left = `${Math.max(0, left)}px`;
    this.tooltip.style.top = `${Math.max(0, y - 10)}px`;
  }

  toBlob(callback) { this.canvas.toBlob(callback, "image/png"); }
}

function drawFrame(ctx, frame, xScale, yScale, xTicks, yTicks, spec) {
  ctx.strokeStyle = THEME.grid;
  ctx.lineWidth = 1;
  ctx.fillStyle = THEME.muted;
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const tick of yTicks) {
    const y = Math.round(yScale(tick.value)) + 0.5;
    ctx.beginPath(); ctx.moveTo(frame.left, y); ctx.lineTo(frame.right, y); ctx.stroke();
    ctx.fillText(tick.label, frame.left - 6, y);
  }
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (const tick of xTicks) {
    const x = Math.round(xScale(tick.value)) + 0.5;
    ctx.beginPath(); ctx.moveTo(x, frame.top); ctx.lineTo(x, frame.bottom); ctx.stroke();
    ctx.fillText(tick.label, x, frame.bottom + 5);
  }
  ctx.strokeStyle = THEME.axis;
  ctx.beginPath(); ctx.moveTo(frame.left, frame.bottom + 0.5); ctx.lineTo(frame.right, frame.bottom + 0.5); ctx.stroke();
  if (spec.yLabel) {
    ctx.save();
    ctx.translate(12, (frame.top + frame.bottom) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(spec.yLabel, 0, 0);
    ctx.restore();
  }
  if (spec.xLabel) {
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillText(spec.xLabel, (frame.left + frame.right) / 2, frame.bottom + 34);
  }
}

/**
 * spec: { series: [{label, color, x, y, dash}], height, xLabel, yLabel,
 *         xTicks: [{value,label}] | null, xFormat, yFormat, zeroLine, yInclude }
 */
export class LineChart extends BaseChart {
  constructor(container, spec) { super(container, spec); this.draw(); }

  draw() {
    const spec = this.spec;
    const { ctx, width, height } = setupCanvas(this.canvas, spec.height ?? 220);
    const frame = { left: 64, right: width - 16, top: 14, bottom: height - (spec.xLabel ? 42 : 26) };
    const xs = spec.series.flatMap((series) => Array.from(series.x));
    const ys = spec.series.flatMap((series) => Array.from(series.y));
    const [xMin, xMax] = spec.xDomain ?? extent(xs);
    let [yMin, yMax] = extent(ys, spec.yInclude ?? []);
    const anchored = spec.yInclude?.includes(yMin);
    if (yMin === yMax) {
      if (anchored) yMax = yMin + 1;
      else { yMin -= Math.abs(yMin || 1) * 0.05; yMax += Math.abs(yMax || 1) * 0.05; }
    }
    const pad = (yMax - yMin) * 0.06;
    yMin -= anchored ? 0 : pad;
    yMax += pad;
    const xScale = (x) => frame.left + ((x - xMin) / (xMax - xMin || 1)) * (frame.right - frame.left);
    const yScale = (y) => frame.bottom - ((y - yMin) / (yMax - yMin || 1)) * (frame.bottom - frame.top);
    const yFormat = spec.yFormat ?? ((value) => formatNumber(value, 3));
    const xFormat = spec.xFormat ?? ((value) => formatNumber(value, 3));
    const yTicks = niceTicks(yMin, yMax, 5).map((value) => ({ value, label: yFormat(value) }));
    const xTicks = spec.xTicks ?? niceTicks(xMin, xMax, 6).map((value) => ({ value, label: xFormat(value) }));
    drawFrame(ctx, frame, xScale, yScale, xTicks.filter((t) => t.value >= xMin && t.value <= xMax), yTicks, spec);
    if (spec.zeroLine && yMin < 0 && yMax > 0) {
      ctx.strokeStyle = THEME.axis;
      ctx.beginPath(); ctx.moveTo(frame.left, yScale(0)); ctx.lineTo(frame.right, yScale(0)); ctx.stroke();
    }
    ctx.save();
    ctx.beginPath(); ctx.rect(frame.left, frame.top - 2, frame.right - frame.left, frame.bottom - frame.top + 4); ctx.clip();
    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    for (const series of spec.series) {
      ctx.strokeStyle = series.color;
      ctx.lineWidth = series.width ?? 2;
      ctx.setLineDash(series.dash ?? []);
      ctx.beginPath();
      let pen = false;
      for (let i = 0; i < series.x.length; i += 1) {
        const y = series.y[i];
        if (!Number.isFinite(y)) { pen = false; continue; }
        const px = xScale(series.x[i]);
        const py = yScale(y);
        if (pen) ctx.lineTo(px, py); else ctx.moveTo(px, py);
        pen = true;
      }
      ctx.stroke();
      if (series.markers) {
        ctx.setLineDash([]);
        for (let i = 0; i < series.x.length; i += 1) {
          if (!Number.isFinite(series.y[i])) continue;
          ctx.beginPath();
          ctx.arc(xScale(series.x[i]), yScale(series.y[i]), 4, 0, Math.PI * 2);
          ctx.fillStyle = series.color;
          ctx.fill();
          ctx.lineWidth = 2;
          ctx.strokeStyle = THEME.surface;
          ctx.stroke();
          ctx.strokeStyle = series.color;
        }
      }
    }
    ctx.restore();
    ctx.setLineDash([]);
    if (!this.hover || this.hover.x < frame.left || this.hover.x > frame.right) { this.tooltip.hidden = true; return; }
    const xValue = xMin + ((this.hover.x - frame.left) / (frame.right - frame.left)) * (xMax - xMin);
    const lines = [];
    let nearestX = null;
    for (const series of spec.series) {
      let best = -1;
      let distance = Infinity;
      for (let i = 0; i < series.x.length; i += 1) {
        const d = Math.abs(series.x[i] - xValue);
        if (d < distance && Number.isFinite(series.y[i])) { distance = d; best = i; }
      }
      if (best < 0) continue;
      nearestX ??= series.x[best];
      const px = xScale(series.x[best]);
      const py = yScale(series.y[best]);
      ctx.beginPath(); ctx.arc(px, py, 4, 0, Math.PI * 2);
      ctx.fillStyle = series.color; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = THEME.surface; ctx.stroke();
      lines.push({ color: series.color, text: `${series.label}: ${yFormat(series.y[best])}${spec.yUnit ? ` ${spec.yUnit}` : ""}` });
    }
    if (nearestX === null) return;
    ctx.strokeStyle = THEME.axis;
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(xScale(nearestX) + 0.5, frame.top); ctx.lineTo(xScale(nearestX) + 0.5, frame.bottom); ctx.stroke();
    this.showTooltip([{ text: (spec.xTooltip ?? xFormat)(nearestX) }, ...lines], this.hover.x, this.hover.y);
  }
}

/**
 * Walter-style climograph: monthly precipitation bars against the right axis
 * and temperature line against the left axis (the conventional 雨温図 layout).
 * spec: { temperature[12], precipitation[12], temperatureRange, precipitationMax, height }
 */
export class ClimographChart extends BaseChart {
  constructor(container, spec) { super(container, spec); this.draw(); }

  draw() {
    const spec = this.spec;
    const { ctx, width, height } = setupCanvas(this.canvas, spec.height ?? 190);
    const frame = { left: 40, right: width - 42, top: 10, bottom: height - 22 };
    // The requested (shared) ranges are widened when the data would fall outside.
    const temperatures = spec.temperature.filter(Number.isFinite);
    const precipitation = spec.precipitation.filter(Number.isFinite);
    const [requestedMin, requestedMax] = spec.temperatureRange;
    const inside = temperatures.every((t) => t >= requestedMin && t <= requestedMax);
    const [tMin, tMax] = inside ? [requestedMin, requestedMax]
      : niceRange(Math.min(requestedMin, ...temperatures), Math.max(requestedMax, ...temperatures), 6);
    const wettest = Math.max(0, ...precipitation);
    const pMax = wettest <= spec.precipitationMax ? spec.precipitationMax : niceRange(0, wettest, 6)[1];
    const slot = (frame.right - frame.left) / 12;
    const xScale = (month) => frame.left + (month - 0.5) * slot;
    const tScale = (t) => frame.bottom - ((t - tMin) / (tMax - tMin)) * (frame.bottom - frame.top);
    const pScale = (p) => frame.bottom - (p / pMax) * (frame.bottom - frame.top);
    const tTicks = niceTicks(tMin, tMax, 5).map((value) => ({ value, label: `${formatNumber(value)}°` }));
    const months = Array.from({ length: 12 }, (_, i) => ({ value: i + 1, label: String(i + 1) }));
    drawFrame(ctx, frame, xScale, tScale, months, tTicks, {});
    ctx.fillStyle = THEME.muted;
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    for (const value of niceTicks(0, pMax, 4)) ctx.fillText(formatNumber(value), frame.right + 6, pScale(value));
    const barWidth = Math.min(24, slot * 0.62);
    ctx.fillStyle = "rgba(57, 135, 229, 0.78)";
    spec.precipitation.forEach((p, i) => {
      if (!Number.isFinite(p)) return;
      const x = xScale(i + 1) - barWidth / 2;
      const top = pScale(Math.max(p, 0));
      if (frame.bottom - top < 0.5) return;
      const radius = Math.max(0, Math.min(4, (frame.bottom - top) / 2, barWidth / 2));
      ctx.beginPath();
      ctx.moveTo(x, frame.bottom);
      ctx.lineTo(x, top + radius);
      ctx.arcTo(x, top, x + radius, top, radius);
      ctx.lineTo(x + barWidth - radius, top);
      ctx.arcTo(x + barWidth, top, x + barWidth, top + radius, radius);
      ctx.lineTo(x + barWidth, frame.bottom);
      ctx.closePath();
      ctx.fill();
    });
    if (tMin < 0 && tMax > 0) {
      ctx.strokeStyle = THEME.axis;
      ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(frame.left, tScale(0)); ctx.lineTo(frame.right, tScale(0)); ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.strokeStyle = "#e66767";
    ctx.lineWidth = 2;
    ctx.lineJoin = "round";
    ctx.beginPath();
    spec.temperature.forEach((t, i) => (i ? ctx.lineTo(xScale(i + 1), tScale(t)) : ctx.moveTo(xScale(i + 1), tScale(t))));
    ctx.stroke();
    spec.temperature.forEach((t, i) => {
      ctx.beginPath(); ctx.arc(xScale(i + 1), tScale(t), 3.5, 0, Math.PI * 2);
      ctx.fillStyle = "#e66767"; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = THEME.surface; ctx.stroke();
    });
    if (!this.hover) { this.tooltip.hidden = true; return; }
    const month = Math.round((this.hover.x - frame.left) / slot + 0.5);
    if (month < 1 || month > 12) { this.tooltip.hidden = true; return; }
    this.showTooltip([
      { text: `${month}月` },
      { color: "#e66767", text: `気温 ${formatNumber(spec.temperature[month - 1], 3)} °C` },
      { color: "#3987e5", text: `降水量 ${formatNumber(spec.precipitation[month - 1], 3)} mm/月` },
    ], this.hover.x, this.hover.y);
  }
}

/** Marching squares on a rectilinear grid; returns segments in data coordinates. */
export function contourSegments(x, y, z, level) {
  const segments = [];
  const at = (i, j) => z[i][j];
  const lerp = (a, b, va, vb) => a + (b - a) * ((level - va) / (vb - va));
  for (let i = 0; i < x.length - 1; i += 1) {
    for (let j = 0; j < y.length - 1; j += 1) {
      const corners = [
        [x[i], y[j], at(i, j)], [x[i + 1], y[j], at(i + 1, j)],
        [x[i + 1], y[j + 1], at(i + 1, j + 1)], [x[i], y[j + 1], at(i, j + 1)],
      ];
      const points = [];
      for (let e = 0; e < 4; e += 1) {
        const [ax, ay, va] = corners[e];
        const [bx, by, vb] = corners[(e + 1) % 4];
        if ((va < level) !== (vb < level)) points.push([lerp(ax, bx, va, vb), lerp(ay, by, va, vb)]);
      }
      if (points.length === 2) segments.push([...points[0], ...points[1]]);
      else if (points.length === 4) segments.push([...points[0], ...points[1]], [...points[2], ...points[3]]);
    }
  }
  return segments;
}

/**
 * Filled contours of z[i][j] over x[i] (ascending) and y[j] (ascending).
 * spec: { x, y, z, levels (ascending, n+1), colors (n), xTicks, yTicks, xLabel,
 *         yLabel, format, tooltip(xValue, yValue, zValue) }
 */
export class ContourChart extends BaseChart {
  constructor(container, spec) { super(container, spec); this.draw(); }

  draw() {
    const spec = this.spec;
    const { ctx, width, height } = setupCanvas(this.canvas, spec.height ?? 300);
    const frame = { left: 64, right: width - 86, top: spec.keyLabel ? 26 : 12, bottom: height - 42 };
    const { x, y, z, levels, colors } = spec;
    const xMin = spec.xDomain?.[0] ?? x[0];
    const xMax = spec.xDomain?.[1] ?? x.at(-1);
    const yMin = y[0];
    const yMax = y.at(-1);
    const xScale = (value) => frame.left + ((value - xMin) / (xMax - xMin)) * (frame.right - frame.left);
    const yScale = (value) => frame.bottom - ((value - yMin) / (yMax - yMin)) * (frame.bottom - frame.top);
    const locate = (values, value) => {
      if (value <= values[0]) return [0, 0];
      if (value >= values.at(-1)) return [values.length - 2, 1];
      let i = 0;
      while (values[i + 1] < value) i += 1;
      return [i, (value - values[i]) / (values[i + 1] - values[i])];
    };
    const sample = (xValue, yValue) => {
      const [i, a] = locate(x, xValue);
      const [j, b] = locate(y, yValue);
      return (1 - a) * (1 - b) * z[i][j] + a * (1 - b) * z[i + 1][j] + (1 - a) * b * z[i][j + 1] + a * b * z[i + 1][j + 1];
    };
    const dpr = this.canvas.width / width;
    const pixelWidth = Math.round((frame.right - frame.left) * dpr);
    const pixelHeight = Math.round((frame.bottom - frame.top) * dpr);
    const image = ctx.createImageData(pixelWidth, pixelHeight);
    const rgb = colors.map((hex) => [1, 3, 5].map((k) => Number.parseInt(hex.slice(k, k + 2), 16)));
    for (let py = 0; py < pixelHeight; py += 1) {
      const yValue = yMax - ((py + 0.5) / pixelHeight) * (yMax - yMin);
      for (let px = 0; px < pixelWidth; px += 1) {
        const xValue = xMin + ((px + 0.5) / pixelWidth) * (xMax - xMin);
        const value = sample(xValue, yValue);
        let band = 0;
        while (band < colors.length - 1 && value > levels[band + 1]) band += 1;
        const offset = (py * pixelWidth + px) * 4;
        image.data[offset] = rgb[band][0];
        image.data[offset + 1] = rgb[band][1];
        image.data[offset + 2] = rgb[band][2];
        image.data[offset + 3] = 255;
      }
    }
    ctx.putImageData(image, Math.round(frame.left * dpr), Math.round(frame.top * dpr));
    const clip = () => { ctx.save(); ctx.beginPath(); ctx.rect(frame.left, frame.top, frame.right - frame.left, frame.bottom - frame.top); ctx.clip(); };
    clip();
    for (const level of levels.slice(1, -1)) {
      ctx.strokeStyle = level === 0 ? "rgba(0,0,0,0.9)" : "rgba(0,0,0,0.3)";
      ctx.lineWidth = level === 0 ? 1.6 : 0.7;
      ctx.beginPath();
      for (const [x0, y0, x1, y1] of contourSegments(x, y, z, level)) {
        ctx.moveTo(xScale(x0), yScale(y0));
        ctx.lineTo(xScale(x1), yScale(y1));
      }
      ctx.stroke();
    }
    ctx.restore();
    ctx.fillStyle = THEME.muted;
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (const tick of spec.yTicks) ctx.fillText(tick.label, frame.left - 6, yScale(tick.value));
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const tick of spec.xTicks) ctx.fillText(tick.label, xScale(tick.value), frame.bottom + 5);
    ctx.strokeStyle = THEME.axis;
    ctx.strokeRect(frame.left + 0.5, frame.top + 0.5, frame.right - frame.left, frame.bottom - frame.top);
    ctx.textBaseline = "bottom";
    if (spec.xLabel) ctx.fillText(spec.xLabel, (frame.left + frame.right) / 2, height - 4);
    if (spec.yLabel) {
      ctx.save(); ctx.translate(12, (frame.top + frame.bottom) / 2); ctx.rotate(-Math.PI / 2);
      ctx.textBaseline = "middle"; ctx.fillText(spec.yLabel, 0, 0); ctx.restore();
    }
    // Colour key.
    const keyLeft = frame.right + 18;
    const bandHeight = (frame.bottom - frame.top) / colors.length;
    colors.forEach((color, band) => {
      ctx.fillStyle = color;
      ctx.fillRect(keyLeft, frame.bottom - (band + 1) * bandHeight, 12, Math.ceil(bandHeight));
    });
    ctx.fillStyle = THEME.muted;
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    const step = Math.max(1, Math.round(levels.length / 7));
    levels.forEach((level, index) => {
      if (index % step && index !== levels.length - 1) return;
      ctx.fillText((spec.format ?? formatNumber)(level), keyLeft + 16, frame.bottom - index * bandHeight);
    });
    if (spec.keyLabel) { ctx.textBaseline = "bottom"; ctx.fillText(spec.keyLabel, keyLeft - 4, frame.top - 8); }
    if (!this.hover || this.hover.x < frame.left || this.hover.x > frame.right
      || this.hover.y < frame.top || this.hover.y > frame.bottom) { this.tooltip.hidden = true; return; }
    const xValue = xMin + ((this.hover.x - frame.left) / (frame.right - frame.left)) * (xMax - xMin);
    const yValue = yMax - ((this.hover.y - frame.top) / (frame.bottom - frame.top)) * (yMax - yMin);
    this.showTooltip(spec.tooltip(xValue, yValue, sample(xValue, yValue)).map((text) => ({ text })), this.hover.x, this.hover.y);
  }
}
