import * as THREE from "/vendor/three.module.js";
import { buildStreamlines, wrapLongitude } from "/streamlines.js";
import * as C from "/climate.js";
import * as F from "/fields.js";
import { ClimographChart, ContourChart, LATITUDE_TICKS, LineChart, SERIES, formatNumber, niceRange } from "/charts.js";

// Colours are handled as plain sRGB values end to end.
THREE.ColorManagement.enabled = false;

const $ = (selector) => document.querySelector(selector);
const ui = {
  run: $("#run-select"), year: $("#year-select"), thresholdRange: $("#threshold-range"),
  thresholdNumber: $("#threshold-number"), datasetState: $("#dataset-state"),
  metaCase: $("#meta-case"), metaMonths: $("#meta-months"), metaLand: $("#meta-land"), metaPoints: $("#meta-points"),
  field: $("#field-select"), levelPicker: $("#level-picker"), level: $("#level-select"),
  coastline: $("#coastline"), sites: $("#sites"), gridPoints: $("#grid-points"),
  streamlines: $("#streamlines"), density: $("#streamline-density"), savePng: $("#save-map"),
  globe: $("#globe"), map: $("#map"), siteLabels: $("#site-labels"), readout: $("#readout"),
  fieldTitle: $("#field-title"), fieldSubtitle: $("#field-subtitle"), legend: $("#legend"),
  transport: $("#transport"), play: $("#play"), timeline: $("#timeline"), frameLabel: $("#frame-label"),
  frameDetail: $("#frame-detail"), rate: $("#rate"), scaleMode: $("#scale-mode"),
  tabs: $("#analysis-tabs"), analysisBody: $("#analysis-body"), analysisState: $("#analysis-state"),
  fatal: $("#fatal"), viewButtons: [...document.querySelectorAll("[data-view-mode]")],
  mapAxes: $("#map-axes"),
};

const CLIMATE_FIELDS = [
  "surface_temperature", "land_temperature", "ocean_temperature", "precipitation", "evaporation",
  "cloud_cover", "surface_pressure", "surface_water", "sea_ice_fraction", "sea_ice_volume",
  "sea_ice_temperature", "snow_fraction",
];
const ANALYSIS_TABS = [
  ["climographs", "雨温図"], ["zonal", "帯状平均"], ["streamfunction", "質量流線関数"],
  ["seaice", "海氷"], ["timeseries", "時系列"], ["comparison", "比較"], ["summary", "要約と CSV"],
];

const state = {
  runs: [], info: null, grid: null, staticFields: {}, daily: null,
  threshold: clampThreshold(Number(localStorage.getItem("bespoke.landThreshold") ?? 0.5)),
  years: [], analysisYear: null, climateInput: null, climate: null, climateFields: null,
  field: null, frames: [], frameIndex: 0, level: 1, levelPressures: [],
  display: null, scaleMode: "auto", scale: null, lockedScale: null, globalScales: new Map(),
  playing: false, rate: 1, lastTick: 0, accumulator: 0, loading: false,
  streamlinesEnabled: true, density: "medium", viewMode: "globe",
  customSite: null, hoverCell: -1, analysisTab: "climographs", referenceRun: null, charts: [],
  token: 0, runToken: 0, climateToken: 0,
};
const cache = new Map();
let cacheBytes = 0;
const views = { globe: null, map: null };

function clampThreshold(value) {
  return Number.isFinite(value) ? Math.min(1, Math.max(0.01, Math.round(value * 100) / 100)) : 0.5;
}

function showFatal(error) {
  console.error(error);
  ui.fatal.textContent = `エラー: ${error.message}`;
  ui.fatal.hidden = false;
}

async function fetchJSON(url) {
  const response = await fetch(url, { cache: "no-store" });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `${response.status} ${response.statusText}`);
  return data;
}

/** Records of float64 values; cached per URL. */
async function fetchRecords(url, recordLength) {
  if (cache.has(url)) {
    const hit = cache.get(url);
    cache.delete(url);
    cache.set(url, hit);
    return hit;
  }
  const response = await fetch(url, { cache: "no-store" });
  if (!response.ok) {
    let message = `${response.status} ${response.statusText}`;
    try { message = (await response.json()).error || message; } catch { /* binary error body */ }
    throw new Error(message);
  }
  const values = new Float64Array(await response.arrayBuffer());
  if (values.length % recordLength) throw new Error(`${url} のデータの長さが格子と合いません`);
  const records = [];
  for (let offset = 0; offset < values.length; offset += recordLength) records.push(values.subarray(offset, offset + recordLength));
  cache.set(url, records);
  cacheBytes += values.byteLength;
  while (cacheBytes > 160e6 && cache.size > 1) {
    const [key, value] = cache.entries().next().value;
    cache.delete(key);
    cacheBytes -= value.reduce((sum, record) => sum + record.byteLength, 0);
  }
  return records;
}

const runUrl = (suffix) => `/api/runs/${encodeURIComponent(state.info.name)}/${suffix}`;

function gridRecords(sampling, field, indices = null, level = null, runName = state.info.name) {
  const params = new URLSearchParams();
  if (indices) params.set("index", indices.join(","));
  if (level !== null) params.set("level", String(level));
  params.set("g", String(state.info.generation));
  const url = `/api/runs/${encodeURIComponent(runName)}/grid/${sampling}/${field}?${params}`;
  return fetchRecords(url, state.grid.pointCount);
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

const DEG = Math.PI / 180;
function globePoint(lonDeg, latDeg, radius = 1) {
  const lon = lonDeg * DEG;
  const lat = latDeg * DEG;
  // The z axis is mirrored so that eastward is screen-right seen from outside.
  return [radius * Math.cos(lat) * Math.cos(lon), radius * Math.sin(lat), -radius * Math.cos(lat) * Math.sin(lon)];
}

/** One flat quad (split into longitude segments on the globe) per native grid cell. */
function buildCellMesh(grid, isMap) {
  const positions = [];
  const sources = [];
  const pushQuad = (cell, lon0, lon1, lat0, lat1) => {
    const corners = isMap
      ? [[lon0 / 180, lat0 / 90, 0], [lon1 / 180, lat0 / 90, 0], [lon1 / 180, lat1 / 90, 0], [lon0 / 180, lat1 / 90, 0]]
      : [globePoint(lon0, lat0, 1.001), globePoint(lon1, lat0, 1.001), globePoint(lon1, lat1, 1.001), globePoint(lon0, lat1, 1.001)];
    for (const corner of [0, 1, 2, 0, 2, 3]) {
      positions.push(...corners[corner]);
      sources.push(cell);
    }
  };
  for (let j = 0; j < grid.nlat; j += 1) {
    const width = 360 / grid.nlon[j];
    const south = grid.latEdges[j];
    const north = grid.latEdges[j + 1];
    for (let k = 0; k < grid.nlon[j]; k += 1) {
      const cell = grid.offsets[j] + k;
      const left = grid.mapLon[cell] - width / 2;
      const right = grid.mapLon[cell] + width / 2;
      if (isMap) {
        if (left < -180) { pushQuad(cell, -180, right, south, north); pushQuad(cell, left + 360, 180, south, north); }
        else if (right > 180) { pushQuad(cell, left, 180, south, north); pushQuad(cell, -180, right - 360, south, north); }
        else pushQuad(cell, left, right, south, north);
      } else {
        const segments = Math.max(1, Math.ceil(width / 3));
        for (let s = 0; s < segments; s += 1) {
          pushQuad(cell, left + (width * s) / segments, left + (width * (s + 1)) / segments, south, north);
        }
      }
    }
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(new Float32Array(positions.length), 3));
  const mesh = new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
  return { mesh, sources: Uint32Array.from(sources) };
}

function lineMaterial(color, opacity, depthTest = true) {
  return new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity, depthTest, depthWrite: false });
}

function addGraticule(scene, isMap) {
  const vertices = [];
  if (isMap) {
    for (let lon = -120; lon <= 120; lon += 60) vertices.push(lon / 180, -1, 0.01, lon / 180, 1, 0.01);
    for (let lat = -60; lat <= 60; lat += 30) vertices.push(-1, lat / 90, 0.01, 1, lat / 90, 0.01);
  } else {
    for (let lat = -60; lat <= 60; lat += 30) {
      for (let lon = 0; lon < 360; lon += 3) vertices.push(...globePoint(lon, lat, 1.002), ...globePoint(lon + 3, lat, 1.002));
    }
    for (let lon = 0; lon < 360; lon += 30) {
      for (let lat = -90; lat < 90; lat += 3) vertices.push(...globePoint(lon, lat, 1.002), ...globePoint(lon, lat + 3, 1.002));
    }
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(vertices, 3));
  const lines = new THREE.LineSegments(geometry, lineMaterial(0xffffff, 0.16));
  lines.renderOrder = 2;
  scene.add(lines);
}

function coastlineObject(segments, isMap) {
  const vertices = [];
  for (const [x0, y0, x1, y1] of segments) {
    if (isMap) { vertices.push(x0 / 180, y0 / 90, 0.02, x1 / 180, y1 / 90, 0.02); continue; }
    const steps = Math.max(1, Math.ceil(Math.max(Math.abs(x1 - x0), Math.abs(y1 - y0)) / 2));
    for (let s = 0; s < steps; s += 1) {
      const a = s / steps;
      const b = (s + 1) / steps;
      vertices.push(...globePoint(x0 + (x1 - x0) * a, y0 + (y1 - y0) * a, 1.003),
        ...globePoint(x0 + (x1 - x0) * b, y0 + (y1 - y0) * b, 1.003));
    }
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(vertices, 3));
  const lines = new THREE.LineSegments(geometry, lineMaterial(0x14181c, 0.95));
  lines.renderOrder = 4;
  return lines;
}

function pointsObject(grid, isMap) {
  const positions = new Float32Array(grid.pointCount * 3);
  for (let i = 0; i < grid.pointCount; i += 1) {
    const p = isMap ? [grid.mapLon[i] / 180, grid.lat[i] / 90, 0.03] : globePoint(grid.lon[i], grid.lat[i], 1.004);
    positions.set(p, i * 3);
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  const points = new THREE.Points(geometry, new THREE.PointsMaterial({
    size: isMap ? 1.6 : 0.008, color: 0x111111, sizeAttenuation: !isMap, transparent: true, opacity: 0.55,
  }));
  points.renderOrder = 5;
  points.visible = ui.gridPoints.checked;
  return points;
}

function createRenderer(container) {
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
  renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 0);
  container.replaceChildren(renderer.domElement);
  return renderer;
}

function setupViews() {
  for (const view of Object.values(views)) {
    if (!view) continue;
    view.observer.disconnect();
    view.renderer.dispose();
  }
  const grid = state.grid;
  for (const key of ["globe", "map"]) {
    const isMap = key === "map";
    const container = ui[key];
    const renderer = createRenderer(container);
    const scene = new THREE.Scene();
    const camera = isMap
      ? new THREE.OrthographicCamera(-1.08, 1.08, 1.08, -1.08, 0.1, 10)
      : new THREE.PerspectiveCamera(36, 1, 0.1, 20);
    camera.position.set(0, 0, isMap ? 2 : 3.3);
    if (!isMap) {
      scene.add(new THREE.Mesh(new THREE.SphereGeometry(0.998, 64, 40), new THREE.MeshBasicMaterial({ color: 0x0a1722 })));
    }
    const cells = buildCellMesh(grid, isMap);
    scene.add(cells.mesh);
    addGraticule(scene, isMap);
    const points = pointsObject(grid, isMap);
    const streamlines = new THREE.Group();
    const overlay = new THREE.Group();
    scene.add(points, streamlines, overlay);
    const view = {
      key, isMap, container, renderer, scene, camera, cells, points, streamlines, overlay, coastline: null,
      tracer: null, yaw: 2.1, pitch: 0.4, distance: 3.3,
    };
    view.resize = () => {
      const width = Math.max(1, container.clientWidth);
      const height = Math.max(1, container.clientHeight);
      renderer.setSize(width, height, false);
      if (!isMap) { camera.aspect = width / height; camera.updateProjectionMatrix(); }
    };
    view.observer = new ResizeObserver(view.resize);
    view.observer.observe(container);
    view.resize();
    installPointer(view);
    views[key] = view;
  }
  setViewMode(state.viewMode);
}

function setViewMode(mode) {
  state.viewMode = mode;
  ui.globe.hidden = mode !== "globe";
  ui.map.hidden = mode !== "map";
  ui.mapAxes.hidden = mode !== "map";
  for (const button of ui.viewButtons) {
    const active = button.dataset.viewMode === mode;
    button.classList.toggle("active", active);
    button.setAttribute("aria-pressed", String(active));
  }
  requestAnimationFrame(() => views[mode]?.resize());
}

/** Longitude/latitude (degrees) under a pointer, or null. */
function pickLonLat(view, event) {
  const rect = view.renderer.domElement.getBoundingClientRect();
  const ndc = new THREE.Vector2(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
  if (view.isMap) {
    const lon = ndc.x * 1.08 * 180;
    const lat = ndc.y * 1.08 * 90;
    return Math.abs(lon) <= 180 && Math.abs(lat) <= 90 ? { lon, lat } : null;
  }
  const raycaster = new THREE.Raycaster();
  raycaster.setFromCamera(ndc, view.camera);
  const hit = raycaster.ray.intersectSphere(new THREE.Sphere(new THREE.Vector3(), 1.001), new THREE.Vector3());
  if (!hit) return null;
  const lat = Math.asin(THREE.MathUtils.clamp(hit.y / 1.001, -1, 1)) / DEG;
  const lon = Math.atan2(-hit.z, hit.x) / DEG;
  return { lon, lat };
}

function installPointer(view) {
  const canvas = view.renderer.domElement;
  let dragging = false;
  let moved = 0;
  let previous = { x: 0, y: 0 };
  canvas.addEventListener("pointerdown", (event) => {
    dragging = !view.isMap;
    moved = 0;
    previous = { x: event.clientX, y: event.clientY };
    canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener("pointermove", (event) => {
    const dx = event.clientX - previous.x;
    const dy = event.clientY - previous.y;
    moved += Math.abs(dx) + Math.abs(dy);
    if (dragging && event.buttons) {
      view.yaw -= dx * 0.006;
      view.pitch = THREE.MathUtils.clamp(view.pitch + dy * 0.006, -1.4, 1.4);
    }
    previous = { x: event.clientX, y: event.clientY };
    const position = pickLonLat(view, event);
    state.hoverCell = position ? C.cellAt(state.grid, position.lon, position.lat) : -1;
    updateReadout(position);
  });
  canvas.addEventListener("pointerup", (event) => {
    dragging = false;
    if (moved > 5) return;
    const position = pickLonLat(view, event);
    if (!position) return;
    const cell = C.cellAt(state.grid, position.lon, position.lat);
    if (cell >= 0) setCustomSite(cell);
  });
  canvas.addEventListener("pointerleave", () => { state.hoverCell = -1; updateReadout(null); });
  if (!view.isMap) {
    canvas.addEventListener("wheel", (event) => {
      event.preventDefault();
      view.distance = THREE.MathUtils.clamp(view.distance + event.deltaY * 0.002, 1.9, 6);
    }, { passive: false });
  }
}

function formatLat(lat) {
  return Math.abs(lat) < 0.05 ? "赤道" : `${lat < 0 ? "南緯" : "北緯"} ${Math.abs(lat).toFixed(1)}°`;
}

function formatLonLat(lon, lat) {
  const east = ((lon + 540) % 360) - 180;
  const longitude = Math.abs(east) < 0.05 || Math.abs(east) > 179.95 ? `経度 ${Math.abs(east).toFixed(1)}°` : `${east < 0 ? "西経" : "東経"} ${Math.abs(east).toFixed(1)}°`;
  return `${formatLat(lat)}・${longitude}`;
}

function updateReadout(position) {
  const cell = state.hoverCell;
  if (!position || cell < 0 || !state.grid) {
    ui.readout.textContent = "ポインターを置くと値を表示・クリックでそのセルの雨温図を追加";
    return;
  }
  const grid = state.grid;
  const parts = [`${formatLonLat(grid.lon[cell], grid.lat[cell])}`, `f_L ${state.staticFields.land_fraction[cell].toFixed(2)}`];
  const field = state.field;
  if (field?.group === "climatology" && state.climate) {
    const spec = state.climateFields[field.id];
    const entry = spec.legend[spec.classify(state.climate, cell)];
    const value = spec.value?.(state.climate, cell);
    const valueText = typeof value === "string" ? value
      : Number.isFinite(value) ? `${spec.valueLabel ? `${spec.valueLabel(state.climate, cell)} ` : ""}${formatNumber(value, 4)}${spec.unit && spec.unit !== "1" ? ` ${spec.unit}` : ""}` : "";
    parts.push(entry ? entry.label : "—");
    if (valueText && !entry?.label.startsWith(valueText)) parts.push(valueText);
  } else if (state.display) {
    const value = state.display.values[cell];
    parts.push(Number.isFinite(value) ? `${formatNumber(value, 5)} ${state.display.unit}` : "タイルなし");
  }
  if (state.climate?.land[cell]) parts.push(`気候型 ${C.KOPPEN_TYPES[state.climate.koppenType[cell]] ?? "—"}`);
  ui.readout.textContent = parts.join("・");
}

// ---------------------------------------------------------------------------
// Colours and overlays
// ---------------------------------------------------------------------------

function applyCellColors(cellColors) {
  for (const view of Object.values(views)) {
    const attribute = view.cells.mesh.geometry.getAttribute("color");
    const { sources } = view.cells;
    for (let vertex = 0; vertex < sources.length; vertex += 1) {
      const offset = sources[vertex] * 3;
      attribute.array[vertex * 3] = cellColors[offset];
      attribute.array[vertex * 3 + 1] = cellColors[offset + 1];
      attribute.array[vertex * 3 + 2] = cellColors[offset + 2];
    }
    attribute.needsUpdate = true;
  }
}

function renderCoastline() {
  const segments = state.climate ? C.coastlineSegments(state.grid, state.climate.land) : [];
  for (const view of Object.values(views)) {
    if (view.coastline) {
      view.overlay.remove(view.coastline);
      view.coastline.geometry.dispose();
    }
    view.coastline = coastlineObject(segments, view.isMap);
    view.coastline.visible = ui.coastline.checked;
    view.overlay.add(view.coastline);
  }
}

function currentSites() {
  if (!state.climate) return [];
  const sites = state.climate.sites.map((site, i) => ({ ...site, label: String(i + 1) }));
  if (state.customSite !== null) sites.unshift({ site: "Selected cell", name: "選択したセル", index: state.customSite, label: "✚", custom: true });
  return sites.filter((site) => site.index >= 0);
}

function renderSiteLabels() {
  const sites = ui.sites.checked ? currentSites() : [];
  ui.siteLabels.replaceChildren(...sites.map((site) => {
    const label = document.createElement("span");
    label.className = site.custom ? "site-label custom" : "site-label";
    label.textContent = site.label;
    label.title = site.name ?? site.site;
    label.dataset.index = String(site.index);
    return label;
  }));
}

function updateSiteLabelPositions() {
  const view = views[state.viewMode];
  if (!view || !state.grid) return;
  const width = view.container.clientWidth;
  const height = view.container.clientHeight;
  const cameraDirection = view.camera.position.clone().normalize();
  for (const label of ui.siteLabels.children) {
    const index = Number(label.dataset.index);
    const point = view.isMap
      ? new THREE.Vector3(state.grid.mapLon[index] / 180, state.grid.lat[index] / 90, 0)
      : new THREE.Vector3(...globePoint(state.grid.lon[index], state.grid.lat[index], 1.01));
    const visible = view.isMap || point.clone().normalize().dot(cameraDirection) > 0.15;
    const projected = point.project(view.camera);
    label.style.transform = `translate(${((projected.x + 1) / 2) * width}px, ${((1 - projected.y) / 2) * height}px) translate(-50%, -50%)`;
    label.hidden = !visible;
  }
}

function renderLegend() {
  const field = state.field;
  if (!field) return;
  if (field.group === "climatology") {
    const spec = state.climateFields[field.id];
    const allowed = spec.legendFilter ? spec.legendFilter(state.climate) : null;
    ui.legend.className = "legend categorical";
    ui.legend.replaceChildren(...spec.legend.filter((entry) => !allowed || allowed.has(entry.label)).map((entry) => {
      const item = document.createElement("span");
      const swatch = document.createElement("i");
      swatch.style.background = entry.color;
      item.append(swatch, document.createTextNode(entry.label));
      return item;
    }));
    return;
  }
  const spec = F.RAW_FIELDS[field.id];
  const scale = state.scale;
  ui.legend.className = "legend continuous";
  const minimum = document.createElement("span");
  minimum.textContent = scale ? `${formatNumber(scale.min, 4)}` : "—";
  const bar = document.createElement("div");
  bar.className = "legend-gradient";
  bar.style.background = F.cssGradient(F.COLORMAPS[spec.kind]);
  const maximum = document.createElement("span");
  maximum.textContent = scale ? `${formatNumber(scale.max, 4)} ${spec.unit}` : "—";
  const missing = document.createElement("span");
  missing.className = "legend-missing";
  const swatch = document.createElement("i");
  swatch.style.background = F.MISSING_COLOR;
  missing.append(swatch, document.createTextNode(spec.mask ? "タイルなし" : ""));
  ui.legend.replaceChildren(minimum, bar, maximum, ...(spec.mask ? [missing] : []));
}

// ---------------------------------------------------------------------------
// Fields
// ---------------------------------------------------------------------------

function populateFieldSelect() {
  const info = state.info;
  const groups = [];
  const option = (group, id, label) => {
    const element = document.createElement("option");
    element.value = `${group}:${id}`;
    element.textContent = label;
    return element;
  };
  const climateGroup = document.createElement("optgroup");
  climateGroup.label = "解析年の気候値";
  for (const id of F.CLIMATOLOGY_ORDER) climateGroup.append(option("climatology", id, state.climateFields[id].label));
  groups.push(climateGroup);
  const monthly = document.createElement("optgroup");
  monthly.label = "月平均";
  for (const id of F.MONTHLY_ORDER) if (info.monthly_fields[id]) monthly.append(option("monthly", id, F.RAW_FIELDS[id].label));
  groups.push(monthly);
  const yearly = document.createElement("optgroup");
  yearly.label = "年次の瞬時値（4月1日）";
  for (const id of F.YEARLY_SURFACE_ORDER) if (info.yearly_surface_fields[id]) yearly.append(option("yearly", id, F.RAW_FIELDS[id].label));
  for (const id of F.YEARLY_LEVEL_ORDER) {
    const available = id === "speed" ? info.yearly_level_fields.u && info.yearly_level_fields.v : info.yearly_level_fields[id];
    if (available) yearly.append(option("yearly", id, `${F.RAW_FIELDS[id].label}（層）`));
  }
  if (yearly.children.length) groups.push(yearly);
  const fixed = document.createElement("optgroup");
  fixed.label = "時間変化しない場";
  for (const id of F.STATIC_ORDER) if (info.static_fields.includes(id)) fixed.append(option("static", id, F.RAW_FIELDS[id].label));
  groups.push(fixed);
  ui.field.replaceChildren(...groups);
}

function isLevelField(field) {
  return field?.group === "yearly" && F.YEARLY_LEVEL_ORDER.includes(field.id);
}

function framesFor(field) {
  const info = state.info;
  if (field.group === "monthly") return info.monthly_fields[field.id];
  if (field.group === "yearly") {
    if (field.id === "speed") return info.yearly_level_fields.u.filter((year) => info.yearly_level_fields.v.includes(year));
    return info.yearly_level_fields[field.id] ?? info.yearly_surface_fields[field.id];
  }
  return [];
}

function frameText(field, index) {
  const frame = state.frames[index];
  const start = C.startMonthOf(state.info.metadata);
  const perYear = Number(state.info.metadata.calendar.months_per_year);
  if (field.group === "monthly") {
    const month = C.calendarMonthOf(frame, start);
    return { label: `${Math.ceil(frame / perYear)} 年目・${month}月`, detail: `月平均・m${String(frame).padStart(4, "0")}` };
  }
  const time = state.info.years.find((year) => year.index === frame)?.time_seconds;
  const yearLength = Number(state.info.metadata.calendar.solar_day_seconds) * Number(state.info.metadata.calendar.days_per_year);
  const year = Number.isFinite(time) ? Math.round(time / yearLength) + 1 : frame;
  return { label: `${year} 年目・${start}月1日`, detail: `瞬時値・y${String(frame).padStart(4, "0")}` };
}

async function setField(value, { keepFrame = false } = {}) {
  const [group, id] = value.split(":");
  const previous = state.field;
  state.field = { group, id };
  ui.field.value = value;
  stopPlayback();
  state.frames = framesFor(state.field);
  if (!keepFrame || previous?.group !== group) {
    if (group === "monthly") {
      const first = state.analysisYear?.months[0];
      state.frameIndex = Math.max(0, state.frames.indexOf(first));
    } else {
      state.frameIndex = Math.max(0, state.frames.length - 1);
    }
  }
  state.frameIndex = Math.min(state.frameIndex, Math.max(0, state.frames.length - 1));
  state.lockedScale = null;
  if (state.scaleMode === "lock") state.scaleMode = "auto";
  ui.scaleMode.value = state.scaleMode;
  const animated = state.frames.length > 0;
  ui.transport.classList.toggle("static", !animated);
  ui.timeline.max = String(Math.max(0, state.frames.length - 1));
  ui.timeline.disabled = !animated;
  ui.play.disabled = !animated;
  ui.scaleMode.disabled = group === "climatology";
  ui.levelPicker.hidden = group !== "yearly";
  ui.streamlines.disabled = group !== "yearly";
  ui.density.disabled = group !== "yearly" || !state.streamlinesEnabled;
  await renderField();
}

async function loadDisplay(field, frameIndex) {
  const spec = F.RAW_FIELDS[field.id];
  const grid = state.grid;
  const frame = state.frames[frameIndex];
  const level = isLevelField(field) ? state.level : null;
  let raw;
  let ice = null;
  if (field.group === "static") {
    [raw] = await gridRecords("static", field.id);
  } else if (field.id === "speed") {
    const [[u], [v]] = await Promise.all([gridRecords("yearly", "u", [frame], level), gridRecords("yearly", "v", [frame], level)]);
    raw = u.map((value, i) => Math.hypot(value, v[i]));
  } else {
    [raw] = await gridRecords(field.group, field.id, [frame], level);
  }
  if (spec.mask === "ice") [ice] = await gridRecords(field.group, "sea_ice_fraction", [frame]);
  return { values: maskAndConvert(spec, raw, ice), unit: spec.unit };
}

function maskAndConvert(spec, raw, ice) {
  const land = state.staticFields.land_fraction;
  const values = new Float64Array(raw.length);
  for (let i = 0; i < raw.length; i += 1) {
    let present = true;
    if (spec.mask === "land") present = land[i] > 0;
    else if (spec.mask === "ocean") present = land[i] < 1;
    else if (spec.mask === "ice") present = land[i] < 1 && ice[i] > 0;
    values[i] = present ? (spec.convert ? spec.convert(raw[i]) : raw[i]) : Number.NaN;
  }
  return values;
}

function percentile(sorted, quantile) {
  if (!sorted.length) return Number.NaN;
  return sorted[Math.min(sorted.length - 1, Math.max(0, Math.ceil(quantile * sorted.length) - 1))];
}

function autoScale(spec, values) {
  const finite = Array.from(values).filter(Number.isFinite).sort((a, b) => a - b);
  if (!finite.length) return { min: 0, max: 1 };
  if (spec.kind === "diverging") {
    const absolute = finite.map(Math.abs).sort((a, b) => a - b);
    const limit = Math.max(percentile(absolute, 0.995), Number.EPSILON);
    return { min: -limit, max: limit };
  }
  const min = spec.zeroBased ? Math.min(0, finite[0]) : percentile(finite, 0.005);
  const max = percentile(finite, 0.995);
  return { min, max: max > min ? max : min + Math.max(Math.abs(min) * 1e-6, Number.EPSILON) };
}

async function globalScale(field, spec) {
  const key = `${state.info.name}:${state.info.generation}:${field.group}:${field.id}:${isLevelField(field) ? state.level : ""}`;
  if (state.globalScales.has(key)) return state.globalScales.get(key);
  let min = Infinity;
  let max = -Infinity;
  let absolute = 0;
  const level = isLevelField(field) ? state.level : null;
  let records;
  if (field.id === "speed") {
    const [u, v] = await Promise.all([gridRecords("yearly", "u", state.frames, level), gridRecords("yearly", "v", state.frames, level)]);
    records = u.map((values, r) => values.map((value, i) => Math.hypot(value, v[r][i])));
  } else {
    records = await gridRecords(field.group, field.id, state.frames, level);
  }
  const ice = spec.mask === "ice" ? await gridRecords(field.group, "sea_ice_fraction", state.frames) : null;
  records.forEach((raw, r) => {
    for (const value of maskAndConvert(spec, raw, ice?.[r])) {
      if (!Number.isFinite(value)) continue;
      min = Math.min(min, value);
      max = Math.max(max, value);
      absolute = Math.max(absolute, Math.abs(value));
    }
  });
  let scale = spec.kind === "diverging" ? { min: -absolute, max: absolute } : { min: spec.zeroBased ? Math.min(0, min) : min, max };
  if (!(scale.max > scale.min)) scale = { min: scale.min, max: scale.min + 1 };
  state.globalScales.set(key, scale);
  return scale;
}

async function renderField() {
  const field = state.field;
  if (!field || !state.grid) return;
  const token = ++state.token;
  const grid = state.grid;
  const cellColors = new Float32Array(grid.pointCount * 3);
  if (field.group === "climatology") {
    if (!state.climate) return;
    const spec = state.climateFields[field.id];
    const colors = spec.legend.map((entry) => F.hexToRgb(entry.color));
    const missing = F.hexToRgb(F.MISSING_COLOR);
    for (let i = 0; i < grid.pointCount; i += 1) cellColors.set(colors[spec.classify(state.climate, i)] ?? missing, i * 3);
    state.display = null;
    state.scale = null;
    ui.fieldTitle.textContent = spec.label;
    ui.fieldSubtitle.textContent = `${state.analysisYear.year} 年目・出力月 ${state.analysisYear.months[0]}–${state.analysisYear.months.at(-1)}・陸 f_L ≥ ${state.threshold.toFixed(2)}`;
    ui.frameLabel.textContent = `${state.analysisYear.year} 年目`;
    ui.frameDetail.textContent = "年平均の気候値";
  } else {
    const spec = F.RAW_FIELDS[field.id];
    state.loading = true;
    let display;
    try {
      display = await loadDisplay(field, state.frameIndex);
    } finally {
      state.loading = false;
    }
    if (token !== state.token) return;
    let scale;
    if (state.scaleMode === "lock" && state.lockedScale) scale = state.lockedScale;
    else if (state.scaleMode === "global" && field.group !== "static") scale = await globalScale(field, spec);
    else scale = autoScale(spec, display.values);
    if (token !== state.token) return;
    state.display = display;
    state.scale = scale;
    const stops = F.COLORMAPS[spec.kind];
    const missing = F.hexToRgb(F.MISSING_COLOR);
    const width = scale.max - scale.min;
    for (let i = 0; i < grid.pointCount; i += 1) {
      const value = display.values[i];
      cellColors.set(Number.isFinite(value) ? F.sampleColormap(stops, (value - scale.min) / width) : missing, i * 3);
    }
    const levelText = isLevelField(field) ? ` · L${String(state.level).padStart(2, "0")} ≈ ${formatPressure(state.levelPressures[state.level - 1])} hPa` : "";
    ui.fieldTitle.textContent = spec.label;
    if (state.frames.length) {
      const text = frameText(field, state.frameIndex);
      ui.frameLabel.textContent = text.label;
      ui.frameDetail.textContent = text.detail;
      ui.fieldSubtitle.textContent = `${text.detail}${levelText}`;
      ui.timeline.value = String(state.frameIndex);
    } else {
      ui.frameLabel.textContent = "時間変化しない場";
      ui.frameDetail.textContent = "";
      ui.fieldSubtitle.textContent = "時間変化しない場";
    }
  }
  applyCellColors(cellColors);
  renderLegend();
  updateReadout(state.hoverCell >= 0 ? {} : null);
  await loadStreamlines(token);
}

function formatPressure(pressureHpa) {
  return pressureHpa >= 100 ? pressureHpa.toFixed(0) : pressureHpa.toFixed(1);
}

// ---------------------------------------------------------------------------
// Streamlines (yearly snapshots)
// ---------------------------------------------------------------------------

function clearStreamlines() {
  for (const view of Object.values(views)) {
    if (!view) continue;
    for (const object of [...view.streamlines.children]) {
      view.streamlines.remove(object);
      object.geometry.dispose();
      object.material.dispose();
    }
    view.tracer = null;
  }
}

function flowPosition(point, isMap) {
  if (isMap) return [wrapLongitude(point.longitude) / Math.PI, point.latitude / (Math.PI / 2), 0.05];
  return globePoint(point.longitude / DEG, point.latitude / DEG, 1.012);
}

async function loadStreamlines(token) {
  clearStreamlines();
  if (state.field?.group !== "yearly" || !state.streamlinesEnabled || !state.frames.length) return;
  const frame = state.frames[state.frameIndex];
  if (!state.info.yearly_level_fields.u?.includes(frame) || !state.info.yearly_level_fields.v?.includes(frame)) return;
  const [[u], [v]] = await Promise.all([
    gridRecords("yearly", "u", [frame], state.level), gridRecords("yearly", "v", [frame], state.level),
  ]);
  if (token !== state.token) return;
  const result = buildStreamlines(state.info.grid, u, v, state.density);
  for (const view of Object.values(views)) {
    const positions = [];
    const colors = [];
    for (const line of result.lines) {
      for (let i = 1; i < line.length; i += 1) {
        const a = line[i - 1];
        const b = line[i];
        if (view.isMap && Math.abs(wrapLongitude(b.longitude) - wrapLongitude(a.longitude)) > Math.PI) continue;
        positions.push(...flowPosition(a, view.isMap), ...flowPosition(b, view.isMap));
        for (const point of [a, b]) {
          const amount = Math.sqrt(Math.min(1, point.speed / Math.max(result.maximumSpeed, Number.EPSILON)));
          colors.push(0.75 + 0.25 * amount, 0.75 + 0.25 * amount, 0.75 + 0.25 * amount);
        }
      }
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
    const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.55, depthWrite: false }));
    lines.renderOrder = 6;
    const tracerGeometry = new THREE.BufferGeometry();
    tracerGeometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array(result.lines.length * 3), 3));
    const tracers = new THREE.Points(tracerGeometry, new THREE.PointsMaterial({
      color: 0xffffff, size: view.isMap ? 3 : 0.02, sizeAttenuation: !view.isMap, depthWrite: false,
    }));
    tracers.renderOrder = 7;
    view.streamlines.add(lines, tracers);
    const speedFactors = result.lines.map((line) => 0.25 + 0.85 * line.reduce((sum, p) => sum + p.speed, 0) / line.length / Math.max(result.maximumSpeed, Number.EPSILON));
    view.tracer = { object: tracers, lines: result.lines, speedFactors };
  }
}

function updateTracers(timestamp) {
  for (const view of Object.values(views)) {
    const tracer = view?.tracer;
    if (!tracer) continue;
    const attribute = tracer.object.geometry.getAttribute("position");
    tracer.lines.forEach((line, index) => {
      const phase = (timestamp * 0.00004 * tracer.speedFactors[index] + index * 0.618034) % 1;
      attribute.array.set(flowPosition(line[Math.min(line.length - 1, Math.floor(phase * line.length))], view.isMap), index * 3);
    });
    attribute.needsUpdate = true;
  }
}

// ---------------------------------------------------------------------------
// Playback
// ---------------------------------------------------------------------------

function stopPlayback() {
  state.playing = false;
  state.accumulator = 0;
  ui.play.classList.remove("playing");
  ui.play.setAttribute("aria-label", "再生");
}

async function stepFrame() {
  if (state.loading || !state.frames.length) return;
  state.frameIndex = (state.frameIndex + 1) % state.frames.length;
  try { await renderField(); } catch (error) { stopPlayback(); showFatal(error); }
}

function animate(timestamp) {
  requestAnimationFrame(animate);
  const view = views[state.viewMode];
  if (view) {
    const globe = views.globe;
    globe.camera.position.set(
      globe.distance * Math.sin(globe.yaw) * Math.cos(globe.pitch),
      globe.distance * Math.sin(globe.pitch),
      globe.distance * Math.cos(globe.yaw) * Math.cos(globe.pitch),
    );
    globe.camera.lookAt(0, 0, 0);
    updateTracers(timestamp);
    view.renderer.render(view.scene, view.camera);
    updateSiteLabelPositions();
  }
  if (state.playing) {
    state.accumulator += (timestamp - state.lastTick) * state.rate;
    if (state.accumulator >= 350) {
      state.accumulator %= 350;
      void stepFrame();
    }
  }
  state.lastTick = timestamp;
}

// ---------------------------------------------------------------------------
// Climate
// ---------------------------------------------------------------------------

function climateOptions(info = state.info, year = state.analysisYear) {
  const metadata = info.metadata;
  const start = C.startMonthOf(metadata);
  return {
    threshold: state.threshold,
    calendarMonths: year.months.map((month) => C.calendarMonthOf(month, start)),
    outputMonths: year.months,
    daysPerMonth: Number(metadata.calendar.days_per_month),
    aHalf: metadata.hybrid_a_half_pa.map(Number),
    bHalf: metadata.hybrid_b_half.map(Number),
    referenceHalfPressure: metadata.reference_half_level_pressure_pa.map(Number),
    gravity: Number(metadata.radiation?.["gravity_acceleration_m_s-2"] ?? 9.80616),
    terrain: info.terrain,
  };
}

async function loadClimateInput() {
  const token = ++state.climateToken;
  const months = state.analysisYear.months;
  ui.analysisState.textContent = `${state.analysisYear.year} 年目を読み込み中…`;
  const nlat = state.grid.nlat;
  const zonalUrl = runUrl(`zonal/zonal_v?index=${months.join(",")}&g=${state.info.generation}`);
  const [records, zonalV] = await Promise.all([
    Promise.all(CLIMATE_FIELDS.map((field) => gridRecords("monthly", field, months))),
    fetchRecords(zonalUrl, nlat * state.info.level_count),
  ]);
  if (token !== state.climateToken) return false;
  state.climateInput = {
    landFraction: state.staticFields.land_fraction,
    surfaceHeight: state.staticFields.surface_height,
    monthly: Object.fromEntries(CLIMATE_FIELDS.map((field, i) => [field, records[i]])),
    zonalV,
  };
  return true;
}

function recomputeClimate() {
  if (!state.climateInput) return;
  state.climate = C.computeClimate(state.grid, state.climateInput, climateOptions());
  let landArea = 0;
  let landCells = 0;
  for (let i = 0; i < state.grid.pointCount; i += 1) {
    if (!state.climate.land[i]) continue;
    landCells += 1;
    landArea += state.grid.weight[i];
  }
  ui.metaLand.textContent = `${landCells.toLocaleString("ja-JP")} セル・${(100 * landArea).toFixed(1)}%`;
  ui.analysisState.textContent = `${state.analysisYear.year} 年目・出力月 ${state.analysisYear.months[0]}–${state.analysisYear.months.at(-1)}・陸 f_L ≥ ${state.threshold.toFixed(2)}`;
  renderCoastline();
  renderSiteLabels();
}

let thresholdFrame = 0;
function setThreshold(value) {
  const threshold = clampThreshold(value);
  ui.thresholdRange.value = String(threshold);
  ui.thresholdNumber.value = threshold.toFixed(2);
  if (threshold === state.threshold && state.climate) return;
  state.threshold = threshold;
  localStorage.setItem("bespoke.landThreshold", String(threshold));
  cancelAnimationFrame(thresholdFrame);
  thresholdFrame = requestAnimationFrame(() => {
    recomputeClimate();
    renderField().catch(showFatal);
    renderAnalysis();
  });
}

function setCustomSite(cell) {
  state.customSite = cell;
  renderSiteLabels();
  if (state.analysisTab !== "climographs") setAnalysisTab("climographs");
  else renderAnalysis();
}

// ---------------------------------------------------------------------------
// Analysis panels
// ---------------------------------------------------------------------------

function destroyCharts() {
  for (const chart of state.charts) chart.destroy();
  state.charts = [];
}

function downloadBlob(blob, name) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function csvCell(value) {
  if (value === null || value === undefined || (typeof value === "number" && !Number.isFinite(value))) return "NA";
  if (typeof value === "number") return String(value);
  const text = String(value);
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function downloadCsv(name, header, rows) {
  const text = [header.join(","), ...rows.map((row) => row.map(csvCell).join(","))].join("\n");
  downloadBlob(new Blob([`${text}\n`], { type: "text/csv" }), `${state.info.name}_${name}`);
}

function button(label, onClick) {
  const element = document.createElement("button");
  element.type = "button";
  element.className = "small-button";
  element.textContent = label;
  element.addEventListener("click", onClick);
  return element;
}

/** A chart card: title, subtitle, legend, body and a PNG button. */
function card(title, subtitle = "", legend = []) {
  const element = document.createElement("article");
  element.className = "chart-card";
  const head = document.createElement("header");
  const heading = document.createElement("div");
  const h = document.createElement("h3");
  h.textContent = title;
  heading.append(h);
  if (subtitle) {
    const p = document.createElement("p");
    p.textContent = subtitle;
    heading.append(p);
  }
  head.append(heading);
  const actions = document.createElement("div");
  actions.className = "card-actions";
  head.append(actions);
  element.append(head);
  if (legend.length) {
    const keys = document.createElement("div");
    keys.className = "series-legend";
    for (const entry of legend) {
      const item = document.createElement("span");
      const key = document.createElement("i");
      key.style.background = entry.color;
      if (entry.dash) key.classList.add("dashed");
      item.append(key, document.createTextNode(entry.label));
      keys.append(item);
    }
    element.append(keys);
  }
  const body = document.createElement("div");
  body.className = "chart-body";
  element.append(body);
  return { element, body, actions, addChart(chart, name) {
    state.charts.push(chart);
    actions.append(button("PNG 保存", () => chart.toBlob((blob) => downloadBlob(blob, `${state.info.name}_${name}.png`))));
    return chart;
  } };
}

function setAnalysisTab(tab) {
  state.analysisTab = tab;
  for (const element of ui.tabs.children) element.classList.toggle("active", element.dataset.tab === tab);
  renderAnalysis();
}

function renderAnalysis() {
  destroyCharts();
  ui.analysisBody.replaceChildren();
  if (!state.climate) return;
  const renderers = {
    climographs: renderClimographs, zonal: renderZonal, streamfunction: renderStreamfunction,
    seaice: renderSeaIce, timeseries: renderTimeSeries, comparison: renderComparison, summary: renderSummary,
  };
  const result = renderers[state.analysisTab]();
  if (result instanceof Promise) result.catch(showFatal);
}

function grid(...children) {
  const element = document.createElement("div");
  element.className = "card-grid";
  element.append(...children);
  return element;
}

function note(text) {
  const element = document.createElement("p");
  element.className = "analysis-note";
  element.textContent = text;
  return element;
}

function renderClimographs() {
  const climate = state.climate;
  const sites = currentSites();
  if (!sites.length) {
    ui.analysisBody.append(note("閾値以上の陸面率を持つセルがありません。"));
    return;
  }
  const records = sites.map((site) => ({ site, record: C.climographOf(climate, site.index) }));
  const presets = records.filter((entry) => !entry.site.custom);
  const pool = presets.length ? presets : records;
  const finite = (values) => values.filter(Number.isFinite);
  const temperatures = finite(pool.flatMap((entry) => entry.record.temperature)).concat(0);
  const precipitation = finite(pool.flatMap((entry) => entry.record.precipitation)).concat(1);
  // Ranges enclose the data of every preset site; the selected cell has its own.
  const sharedRange = niceRange(Math.min(...temperatures), Math.max(...temperatures), 6);
  const sharedMax = niceRange(0, Math.max(...precipitation), 6)[1];
  const cards = records.map(({ site, record }) => {
    const grid = state.grid;
    const height = state.staticFields.surface_height[site.index];
    const label = site.custom ? "✚ 選択したセル" : `${site.label}. ${site.name ?? site.site}`;
    const details = [
      `${formatLonLat(grid.lon[site.index], grid.lat[site.index])}・${Math.round(height)} m・f_L ${state.staticFields.land_fraction[site.index].toFixed(2)}`,
      `${formatNumber(record.annualTemperature, 3)} °C・${Math.round(record.annualPrecipitation)} mm/年`,
    ];
    if (record.temperatureSource === "surface") details.push("陸タイルがないので格子平均の地表温度");
    if (site.earth) details.push(`地球での観測 ${site.earth}`);
    if (Number.isFinite(site.snapDistance)) details.push(`吸着距離 ${site.snapDistance.toFixed(1)}°`);
    const isLand = climate.land[site.index];
    const c = card(`${label} — ${isLand ? record.koppenType : "陸ではない"}`, details.join("・"));
    const custom = site.custom;
    const own = finite(record.temperature);
    const temperatureRange = custom ? niceRange(Math.min(0, ...own), Math.max(0, ...own), 6) : sharedRange;
    const precipitationMax = custom ? niceRange(0, Math.max(1, ...finite(record.precipitation)), 6)[1] : sharedMax;
    c.addChart(new ClimographChart(c.body, {
      temperature: record.temperature, precipitation: record.precipitation, temperatureRange, precipitationMax,
    }), site.custom ? "climograph_selected" : `climograph_${site.label}`);
    if (custom) c.actions.prepend(button("削除", () => { state.customSite = null; renderSiteLabels(); renderAnalysis(); }));
    return c.element;
  });
  ui.analysisBody.append(
    note("1〜12月の月平均の気温（折れ線、左軸、°C）と降水量（棒、右軸、mm/月）。気温は陸温度で、陸タイルのないセルだけ格子平均の地表温度を使う。代表地点は両軸を共有し、陸面率が閾値以上のセルのうち最も近いセルに吸着する。地図をクリックすると任意のセルを加えられる。"),
    grid(...cards),
  );
}

function renderZonal() {
  const zonal = state.climate.zonal;
  const latitude = Array.from(zonal.latitude);
  const xSpec = { xTicks: LATITUDE_TICKS, xDomain: [-90, 90], xTooltip: (x) => formatLat(x) };
  const t = card("年平均の地表温度", "陸と海の両方を含む帯状平均（°C）");
  t.addChart(new LineChart(t.body, { ...xSpec, yUnit: "°C", zeroLine: true, series: [
    { label: "地表温度", color: SERIES[0], x: latitude, y: Array.from(zonal.surfaceTemperatureK, (v) => v - 273.15) },
  ] }), "zonal_temperature");
  const p = card("年平均の降水量と蒸発量", "mm/日", [
    { color: SERIES[0], label: "降水量" }, { color: SERIES[1], label: "蒸発量" }]);
  p.addChart(new LineChart(p.body, { ...xSpec, yUnit: "mm/日", yInclude: [0], series: [
    { label: "降水量", color: SERIES[0], x: latitude, y: Array.from(zonal.precipitation) },
    { label: "蒸発量", color: SERIES[1], x: latitude, y: Array.from(zonal.evaporation) },
  ] }), "zonal_water");
  const c = card("年平均の実効雲量", "診断した鉛直積算の雲量");
  c.addChart(new LineChart(c.body, { ...xSpec, yInclude: [0], series: [
    { label: "雲量", color: SERIES[0], x: latitude, y: Array.from(zonal.cloudCover) },
  ] }), "zonal_cloud");
  ui.analysisBody.append(grid(t.element, p.element, c.element));
}

const STREAM_COLORS = ["#002F70", "#014287", "#2056A2", "#3369BF", "#567DCA", "#7390D4", "#8DA4DD", "#A5B6E6",
  "#BBC8ED", "#CFD8F3", "#E1E7F7", "#F1F3F8", "#F9F1F1", "#F9E1E1", "#F6CFCF", "#F0BBBB", "#E7A6A6", "#DD9090",
  "#D07979", "#C26161", "#B14949", "#953838", "#792727", "#5F1415"];

function renderStreamfunction() {
  const climate = state.climate;
  const pressures = climate.referencePressureHpa;
  const nlev = pressures.length;
  const order = Array.from({ length: nlev }, (_, k) => nlev - 1 - k); // bottom to top
  const y = order.map((k) => -Math.log10(pressures[k]));
  const z = climate.streamfunction.psi.map((row) => order.map((k) => row[k] / 1e9));
  const maximum = Math.max(...z.flat().map(Math.abs));
  const limit = Math.max(1, Math.ceil(maximum / 10) * 10);
  const levels = Array.from({ length: 25 }, (_, i) => -limit + (2 * limit * i) / 24);
  const ticks = [1000, 850, 700, 500, 300, 200, 100, 50, 10, 2]
    .filter((p) => p >= Math.min(...pressures) && p <= Math.max(...pressures))
    .map((p) => ({ value: -Math.log10(p), label: String(p) }));
  const c = card("年平均の子午面質量流線関数",
    `10⁹ kg/s・月平均の帯状平均 v と p_s から計算・縦軸は層の基準気圧・最大 ${formatNumber(Math.max(...z.flat()), 3)}、最小 ${formatNumber(Math.min(...z.flat()), 3)}`);
  c.addChart(new ContourChart(c.body, {
    x: Array.from(climate.zonal.latitude), y, z, levels, colors: STREAM_COLORS, height: 340,
    xDomain: [-90, 90], xTicks: LATITUDE_TICKS, yTicks: ticks, xLabel: "緯度", yLabel: "気圧（hPa）",
    keyLabel: "10⁹ kg/s", format: (v) => formatNumber(v, 3),
    tooltip: (lat, yValue, value) => [formatLat(lat), `${formatPressure(10 ** -yValue)} hPa`, `ψ ${formatNumber(value, 3)} × 10⁹ kg/s`],
  }), "mass_streamfunction");
  c.actions.prepend(button("CSV 保存", downloadStreamfunctionCsv));
  ui.analysisBody.append(c.element);
}

function renderSeaIce() {
  const rows = state.climate.seaIceMonthly;
  const months = Array.from({ length: 12 }, (_, i) => i + 1);
  const series = (region, key, scale) => months.map((month) => {
    const row = rows.find((entry) => entry.region === region && entry.calendarMonth === month);
    return row ? row[key] / scale : Number.NaN;
  });
  const xSpec = {
    xDomain: [0.6, 12.4], xTicks: months.map((m) => ({ value: m, label: `${m}月` })),
    xTooltip: (m) => `${Math.round(m)}月`,
  };
  const legend = [{ color: SERIES[0], label: "北半球" }, { color: SERIES[1], label: "南半球" }];
  const make = (title, key, scale, unit, name) => {
    const c = card(title, unit, legend);
    c.addChart(new LineChart(c.body, { ...xSpec, yInclude: [0], yUnit: unit, series: [
      { label: "北半球", color: SERIES[0], x: months, y: series("north", key, scale), markers: true },
      { label: "南半球", color: SERIES[1], x: months, y: series("south", key, scale), markers: true },
    ] }), name);
    return c;
  };
  const area = make("海氷面積", "ice_area_m2", 1e12, "10¹² m²", "sea_ice_area");
  const volume = make("海氷体積", "ice_volume_m3", 1e13, "10¹³ m³", "sea_ice_volume");
  area.actions.prepend(button("CSV 保存", downloadSeaIceCsv));
  const noIce = rows.every((row) => !(row.areaM2 > 0));
  ui.analysisBody.append(
    ...(noIce ? [note(`${state.analysisYear.year} 年目には海氷がありません。`)] : []),
    note("月平均の海氷面積率 A と体積 V を、元の格子セルの面積に海の割合を掛けた重みで積分したもの。地図は「表示する場」の「解析年の気候値」から海氷・積雪の項目を選ぶ。"),
    grid(area.element, volume.element),
  );
}

function dailySeries(daily, column) {
  const values = daily.data[column];
  return values ? values.map((value) => (value === null ? Number.NaN : value)) : null;
}

function renderTimeSeries() {
  const daily = state.daily;
  const perYear = Number(state.info.metadata.calendar.days_per_year);
  const years = daily.data.simulation_day.map((day) => day / perYear);
  const xSpec = { xLabel: "経過年数", xTooltip: (x) => `経過 ${x.toFixed(2)} 年` };
  const col = (name) => dailySeries(daily, name);
  const toa = years.map((_, i) => col("incoming_shortwave_w_m-2")[i] - col("reflected_shortwave_w_m-2")[i] - col("outgoing_longwave_w_m-2")[i]);
  const cards = [];
  const make = (title, unit, series, extra, name) => {
    const c = card(title, unit, series.length > 1 ? series.map((s) => ({ color: s.color, label: s.label })) : []);
    c.addChart(new LineChart(c.body, { ...xSpec, yUnit: unit, ...extra, series: series.map((s) => ({ ...s, x: years, width: 1.5 })) }), name);
    cards.push(c.element);
  };
  make("全球平均の地表温度", "K", [
    { label: "全体", color: SERIES[0], y: col("mean_surface_temperature_k") },
    { label: "陸", color: SERIES[1], y: col("mean_land_surface_temperature_k") },
    { label: "海", color: SERIES[2], y: col("mean_ocean_surface_temperature_k") },
  ], {}, "global_surface_temperature");
  make("全球平均の降水量と蒸発量", "mm/日", [
    { label: "降水量", color: SERIES[0], y: col("precipitation_mm_day-1") },
    { label: "蒸発量", color: SERIES[1], y: col("evaporation_mm_day-1") },
  ], {}, "global_precipitation_evaporation");
  make("可降水量", "kg/m²", [{ label: "可降水量", color: SERIES[0], y: col("precipitable_water_kg_m-2") }], {}, "precipitable_water");
  make("大気上端の正味放射", "W/m²", [{ label: "正味放射", color: SERIES[0], y: toa }], { zeroLine: true }, "toa_net_radiation");
  make("全球の海氷面積", "10¹² m²", [{ label: "面積", color: SERIES[0], y: col("sea_ice_area_m2").map((v) => v / 1e12) }], { yInclude: [0] }, "sea_ice_area_daily");
  make("全球の海氷体積", "10¹³ m³", [{ label: "体積", color: SERIES[0], y: col("sea_ice_volume_m3").map((v) => v / 1e13) }], { yInclude: [0] }, "sea_ice_volume_daily");
  make("海氷の厚さ（氷面積で重み付け）", "m", [{ label: "厚さ", color: SERIES[0], y: col("mean_sea_ice_thickness_m") }], { yInclude: [0] }, "sea_ice_thickness_daily");
  make("陸の積雪（水当量）", "kg/m²", [{ label: "積雪", color: SERIES[0], y: col("land_snow_water_kg_m-2") }], { yInclude: [0] }, "land_snow_daily");
  ui.analysisBody.append(note("daily_global.csv の全期間の日平均の全球平均（スピンアップの確認）。"), grid(...cards));
}

async function renderComparison() {
  const others = state.runs.filter((run) => run.name !== state.info.name);
  if (!others.length) {
    ui.analysisBody.append(note("output/ に比較できる別の陸海ケースがありません。"));
    return;
  }
  if (!others.some((run) => run.name === state.referenceRun)) {
    state.referenceRun = (others.find((run) => run.terrain !== state.info.terrain && run.truncation === state.info.truncation)
      ?? others.find((run) => run.terrain !== state.info.terrain) ?? others[0]).name;
  }
  const picker = document.createElement("label");
  picker.className = "inline-picker";
  picker.append("比較するケース ");
  const select = document.createElement("select");
  select.replaceChildren(...others.map((run) => {
    const option = document.createElement("option");
    option.value = run.name;
    option.textContent = `${run.case_name}（${run.terrain === "earth" ? "地球の地形" : "解析的な地形"}、T${run.truncation}）`;
    return option;
  }));
  select.value = state.referenceRun;
  select.addEventListener("change", () => { state.referenceRun = select.value; renderAnalysis(); });
  picker.append(select);
  const status = note("比較するケースを読み込み中…");
  ui.analysisBody.append(picker, status);
  const tab = state.analysisTab;
  const name = state.referenceRun;
  const reference = await loadReference(name);
  if (state.analysisTab !== tab || state.referenceRun !== name) return;
  status.textContent = `色つき: ${state.info.case_name} の ${state.analysisYear.year} 年目。灰色（帯状平均は破線）: ${name} の ${reference.year.year} 年目。`;
  const legend = [{ color: SERIES[0], label: state.info.case_name }, { color: "#8b9aa6", label: name, dash: true }];
  const own = state.climate.zonal;
  const latitude = Array.from(own.latitude);
  const xSpec = { xTicks: LATITUDE_TICKS, xDomain: [-90, 90], xTooltip: (x) => formatLat(x) };
  const t = card("帯状平均の地表温度", "°C", legend);
  t.addChart(new LineChart(t.body, { ...xSpec, yUnit: "°C", series: [
    { label: state.info.case_name, color: SERIES[0], x: latitude, y: Array.from(own.surfaceTemperatureK, (v) => v - 273.15) },
    { label: name, color: "#8b9aa6", dash: [5, 4], x: reference.latitude, y: reference.temperature },
  ] }), "comparison_zonal_temperature");
  const p = card("帯状平均の降水量", "mm/日", legend);
  p.addChart(new LineChart(p.body, { ...xSpec, yUnit: "mm/日", yInclude: [0], series: [
    { label: state.info.case_name, color: SERIES[0], x: latitude, y: Array.from(own.precipitation) },
    { label: name, color: "#8b9aa6", dash: [5, 4], x: reference.latitude, y: reference.precipitation },
  ] }), "comparison_zonal_precipitation");
  const perYear = Number(state.info.metadata.calendar.days_per_year);
  const daily = (source, key) => source.data[key].map((v) => (v === null ? Number.NaN : v));
  const toa = (source) => source.data.simulation_day.map((_, i) => daily(source, "incoming_shortwave_w_m-2")[i]
    - daily(source, "reflected_shortwave_w_m-2")[i] - daily(source, "outgoing_longwave_w_m-2")[i]);
  const years = (source) => source.data.simulation_day.map((day) => day / perYear);
  const xDaily = { xLabel: "経過年数", xTooltip: (x) => `経過 ${x.toFixed(2)} 年` };
  const g = card("全球平均の地表温度", "K", legend);
  g.addChart(new LineChart(g.body, { ...xDaily, yUnit: "K", series: [
    { label: state.info.case_name, color: SERIES[0], x: years(state.daily), y: daily(state.daily, "mean_surface_temperature_k"), width: 1.5 },
    { label: name, color: "#8b9aa6", x: years(reference.daily), y: daily(reference.daily, "mean_surface_temperature_k"), width: 1.2 },
  ] }), "comparison_global_temperature");
  const r = card("大気上端の正味放射", "W/m²", legend);
  r.addChart(new LineChart(r.body, { ...xDaily, yUnit: "W/m²", zeroLine: true, series: [
    { label: state.info.case_name, color: SERIES[0], x: years(state.daily), y: toa(state.daily), width: 1.3 },
    { label: name, color: "#8b9aa6", x: years(reference.daily), y: toa(reference.daily), width: 1.1 },
  ] }), "comparison_toa");
  ui.analysisBody.append(grid(t.element, p.element, g.element, r.element));
}

const referenceCache = new Map();
async function loadReference(name) {
  if (referenceCache.has(name)) return referenceCache.get(name);
  const base = `/api/runs/${encodeURIComponent(name)}`;
  const [info, daily] = await Promise.all([fetchJSON(`${base}/metadata`), fetchJSON(`${base}/daily`)]);
  const referenceGrid = C.createGrid(info);
  const years = C.completeYears(info.months, Number(info.metadata.calendar.months_per_year));
  const year = years.find((entry) => entry.year === state.analysisYear.year) ?? years.at(-1);
  const load = (field) => fetchRecords(`${base}/grid/monthly/${field}?index=${year.months.join(",")}&g=${info.generation}`, referenceGrid.pointCount);
  const [temperature, precipitation] = await Promise.all([load("surface_temperature"), load("precipitation")]);
  const result = {
    year, daily,
    latitude: Array.from(referenceGrid.latitude),
    temperature: Array.from(C.zonalAnnualMean(referenceGrid, temperature), (v) => v - 273.15),
    precipitation: Array.from(C.zonalAnnualMean(referenceGrid, precipitation)),
  };
  referenceCache.set(name, result);
  return result;
}

function renderSummary() {
  const climate = state.climate;
  const table = document.createElement("table");
  table.className = "summary-table";
  const head = document.createElement("thead");
  head.innerHTML = "<tr><th>項目</th><th>値</th><th>CSV のキー</th></tr>";
  const body = document.createElement("tbody");
  for (const row of climate.summary) {
    const tr = document.createElement("tr");
    for (const text of [row.label, formatNumber(row.value, 5), row.metric]) {
      const td = document.createElement("td");
      td.textContent = text;
      tr.append(td);
    }
    body.append(tr);
  }
  table.append(head, body);
  const downloads = document.createElement("div");
  downloads.className = "download-row";
  downloads.append(
    button("要約 CSV", () => downloadCsv("analysis_summary.csv", ["metric", "value"], climate.summary.map((row) => [row.metric, row.value]))),
    button("全セルの年平均と気候区分 CSV", downloadCellsCsv),
    button("代表地点 CSV", downloadSitesCsv),
    button("帯状平均 CSV", downloadZonalCsv),
    button("質量流線関数 CSV", downloadStreamfunctionCsv),
    button("海氷の月別集計 CSV", downloadSeaIceCsv),
  );
  const sites = document.createElement("table");
  sites.className = "summary-table";
  sites.innerHTML = "<thead><tr><th>#</th><th>地点</th><th>セル</th><th>吸着距離</th><th>気候群</th><th>気候型</th><th>地球での観測</th><th>気温（°C）</th><th>降水量（mm/年）</th></tr></thead>";
  const siteBody = document.createElement("tbody");
  climate.sites.forEach((site, n) => {
    const tr = document.createElement("tr");
    const i = site.index;
    const cells = i < 0 ? [n + 1, site.name ?? site.site, "—", "—", "—", "—", site.earth ?? "", "—", "—"] : [
      n + 1, site.name ?? site.site, formatLonLat(state.grid.lon[i], state.grid.lat[i]), `${site.snapDistance.toFixed(1)}°`,
      C.KOPPEN_GROUPS[climate.koppenGroup[i]], C.KOPPEN_TYPES[climate.koppenType[i]], site.earth ?? "",
      formatNumber(climate.landTemperatureC[i], 3), Math.round(climate.precipitation[i]),
    ];
    for (const text of cells) {
      const td = document.createElement("td");
      td.textContent = String(text);
      tr.append(td);
    }
    siteBody.append(tr);
  });
  sites.append(siteBody);
  ui.analysisBody.append(downloads, grid(wrapTable("解析年の要約", table), wrapTable("代表地点", sites)));
}

function wrapTable(title, table) {
  const c = card(title);
  c.body.classList.add("table-body");
  c.body.append(table);
  return c.element;
}

function downloadCellsCsv() {
  const c = state.climate;
  const g = state.grid;
  const f = state.staticFields.land_fraction;
  const rows = [];
  for (let i = 0; i < g.pointCount; i += 1) {
    rows.push([
      g.lon[i], g.lat[i], f[i], state.staticFields.surface_height[i], c.land[i] ? 1 : 0,
      c.surfaceTemperatureC[i], f[i] > 0 ? c.landTemperatureC[i] : null, f[i] < 1 ? c.oceanTemperatureC[i] : null,
      c.precipitation[i], c.evaporation[i], c.cloudCover[i], c.surfaceWater[i],
      c.land[i] ? C.KOPPEN_GROUPS[c.koppenGroup[i]] : null, c.land[i] ? C.KOPPEN_TYPES[c.koppenType[i]] : null,
    ]);
  }
  downloadCsv("final_year_koppen_groups.csv", [
    "longitude_deg", "latitude_deg", "land_fraction", "surface_height_m", "is_land",
    "annual_mean_surface_temperature_c", "annual_mean_land_temperature_c", "annual_mean_ocean_temperature_c",
    "annual_precipitation_mm", "annual_evaporation_mm", "annual_mean_cloud_cover", "annual_mean_surface_water_kg_m2",
    "koppen_group", "koppen_type",
  ], rows);
}

function downloadSitesCsv() {
  const c = state.climate;
  const g = state.grid;
  downloadCsv("representative_land_sites.csv", [
    "site", "target_longitude_deg", "target_latitude_deg", "longitude_deg", "latitude_deg", "snap_distance_deg",
    "land_fraction", "surface_height_m", "koppen_group", "koppen_type", "earth_reference_koppen",
    "annual_mean_temperature_c", "annual_precipitation_mm", "annual_mean_surface_water_kg_m2",
  ], c.sites.map((site) => {
    const i = site.index;
    if (i < 0) return [site.site, site.lon, site.lat];
    return [site.site, site.lon, site.lat, g.lon[i], g.lat[i], site.snapDistance, state.staticFields.land_fraction[i],
      state.staticFields.surface_height[i], C.KOPPEN_GROUPS[c.koppenGroup[i]], C.KOPPEN_TYPES[c.koppenType[i]],
      site.earth ?? null, c.landTemperatureC[i], c.precipitation[i], c.surfaceWater[i]];
  }));
}

function downloadZonalCsv() {
  const z = state.climate.zonal;
  downloadCsv("final_year_zonal_means.csv", [
    "latitude_deg", "annual_mean_surface_temperature_k", "annual_mean_precipitation_mm_day",
    "annual_mean_evaporation_mm_day", "annual_mean_cloud_cover",
  ], Array.from(z.latitude, (lat, j) => [lat, z.surfaceTemperatureK[j], z.precipitation[j], z.evaporation[j], z.cloudCover[j]]));
}

function downloadStreamfunctionCsv() {
  const c = state.climate;
  const rows = [];
  c.referencePressureHpa.forEach((pressure, k) => {
    c.zonal.latitude.forEach((lat, j) => rows.push([lat, k + 1, pressure, c.streamfunction.psi[j][k]]));
  });
  downloadCsv("final_year_mass_streamfunction.csv", ["latitude_deg", "level", "reference_pressure_hpa", "mass_streamfunction_kg_s"], rows);
}

function downloadSeaIceCsv() {
  downloadCsv("final_year_sea_ice_monthly.csv",
    ["region", "output_month", "calendar_month", "month", "ice_area_m2", "ice_volume_m3", "mean_thickness_m"],
    state.climate.seaIceMonthly.map((row) => [row.region, row.outputMonth, row.calendarMonth,
      C.MONTH_ABBREVIATIONS[row.calendarMonth - 1], row.areaM2, row.volumeM3, row.thicknessM]));
}

// ---------------------------------------------------------------------------
// Run selection
// ---------------------------------------------------------------------------

async function selectRun(name) {
  const token = ++state.runToken;
  stopPlayback();
  ui.datasetState.textContent = "読み込み中";
  const base = `/api/runs/${encodeURIComponent(name)}`;
  const [info, daily] = await Promise.all([fetchJSON(`${base}/metadata`), fetchJSON(`${base}/daily`)]);
  if (token !== state.runToken) return;
  state.info = info;
  state.daily = daily;
  ui.run.value = name;
  state.grid = C.createGrid(info);
  state.globalScales.clear();
  state.customSite = null;
  state.climate = null;
  state.climateInput = null;
  const statics = await Promise.all(info.static_fields.map((field) => gridRecords("static", field)));
  if (token !== state.runToken) return;
  state.staticFields = Object.fromEntries(info.static_fields.map((field, i) => [field, statics[i][0]]));
  state.climateFields = F.climatologyFields({ maskingWaterEquivalent: Number(info.metadata.snow?.["masking_water_equivalent_kg_m-2"] ?? 50) });
  const half = info.metadata.reference_half_level_pressure_pa.map(Number);
  state.levelPressures = half.slice(0, -1).map((value, k) => 0.5 * (value + half[k + 1]) / 100);
  state.level = state.levelPressures.reduce((best, p, k) => (Math.abs(Math.log(p / 500)) < Math.abs(Math.log(state.levelPressures[best - 1] / 500)) ? k + 1 : best), 1);
  ui.level.replaceChildren(...state.levelPressures.map((pressure, k) => {
    const option = document.createElement("option");
    option.value = String(k + 1);
    option.textContent = `L${String(k + 1).padStart(2, "0")} · ≈ ${formatPressure(pressure)} hPa`;
    return option;
  }));
  ui.level.value = String(state.level);
  state.years = C.completeYears(info.months, Number(info.metadata.calendar.months_per_year));
  state.analysisYear = state.years.at(-1);
  ui.year.replaceChildren(...state.years.map((year) => {
    const option = document.createElement("option");
    option.value = String(year.year);
    option.textContent = `${year.year} 年目（出力月 ${year.months[0]}–${year.months.at(-1)}）`;
    return option;
  }));
  ui.year.value = String(state.analysisYear.year);
  ui.metaCase.textContent = `${info.terrain === "earth" ? "地球の地形" : "解析的な地形"}・T${info.truncation}`;
  ui.metaMonths.textContent = `月平均 ${info.months.length} か月・瞬時値 ${info.years.length} 回`;
  ui.metaPoints.textContent = state.grid.pointCount.toLocaleString("ja-JP");
  setupViews();
  populateFieldSelect();
  if (!(await loadClimateInput()) || token !== state.runToken) return;
  recomputeClimate();
  const keep = state.field && [...ui.field.options].some((option) => option.value === `${state.field.group}:${state.field.id}`);
  await setField(keep ? `${state.field.group}:${state.field.id}` : "climatology:koppen_type");
  renderAnalysis();
  ui.datasetState.textContent = "準備完了";
}

async function selectYear(yearNumber) {
  state.analysisYear = state.years.find((year) => year.year === yearNumber) ?? state.analysisYear;
  if (!(await loadClimateInput())) return;
  recomputeClimate();
  await renderField();
  renderAnalysis();
}

async function init() {
  ui.thresholdRange.value = String(state.threshold);
  ui.thresholdNumber.value = state.threshold.toFixed(2);
  ui.tabs.replaceChildren(...ANALYSIS_TABS.map(([id, label]) => {
    const element = document.createElement("button");
    element.type = "button";
    element.dataset.tab = id;
    element.textContent = label;
    element.classList.toggle("active", id === state.analysisTab);
    element.addEventListener("click", () => setAnalysisTab(id));
    return element;
  }));
  const data = await fetchJSON("/api/runs");
  state.runs = data.runs;
  if (!state.runs.length) throw new Error("output/ に陸海ケース（海氷・雪を含む最新形式）の計算結果がありません");
  ui.run.replaceChildren(...state.runs.map((run) => {
    const option = document.createElement("option");
    option.value = run.name;
    option.textContent = run.case_name;
    return option;
  }));
  ui.run.addEventListener("change", () => selectRun(ui.run.value).catch(showFatal));
  ui.year.addEventListener("change", () => selectYear(Number(ui.year.value)).catch(showFatal));
  ui.thresholdRange.addEventListener("input", () => setThreshold(Number(ui.thresholdRange.value)));
  ui.thresholdNumber.addEventListener("change", () => setThreshold(Number(ui.thresholdNumber.value)));
  ui.field.addEventListener("change", () => setField(ui.field.value, { keepFrame: true }).catch(showFatal));
  ui.level.addEventListener("change", () => { state.level = Number(ui.level.value); renderField().catch(showFatal); });
  ui.coastline.addEventListener("change", () => { for (const view of Object.values(views)) if (view?.coastline) view.coastline.visible = ui.coastline.checked; });
  ui.sites.addEventListener("change", renderSiteLabels);
  ui.gridPoints.addEventListener("change", () => { for (const view of Object.values(views)) if (view) view.points.visible = ui.gridPoints.checked; });
  ui.streamlines.addEventListener("change", () => {
    state.streamlinesEnabled = ui.streamlines.checked;
    ui.density.disabled = !state.streamlinesEnabled || state.field?.group !== "yearly";
    loadStreamlines(state.token).catch(showFatal);
  });
  ui.density.addEventListener("change", () => { state.density = ui.density.value; loadStreamlines(state.token).catch(showFatal); });
  ui.savePng.addEventListener("click", () => {
    const view = views[state.viewMode];
    view.renderer.render(view.scene, view.camera);
    view.renderer.domElement.toBlob((blob) => downloadBlob(blob, `${state.info.name}_${state.field.id}_${state.viewMode}.png`));
  });
  ui.play.addEventListener("click", () => {
    if (state.playing) { stopPlayback(); return; }
    state.playing = true;
    state.accumulator = 0;
    ui.play.classList.add("playing");
    ui.play.setAttribute("aria-label", "一時停止");
  });
  ui.timeline.addEventListener("input", () => {
    stopPlayback();
    state.frameIndex = Number(ui.timeline.value);
    renderField().catch(showFatal);
  });
  ui.rate.addEventListener("change", () => { state.rate = Number(ui.rate.value); });
  ui.scaleMode.addEventListener("change", () => {
    state.scaleMode = ui.scaleMode.value;
    state.lockedScale = state.scaleMode === "lock" ? state.scale : null;
    renderField().catch(showFatal);
  });
  for (const element of ui.viewButtons) element.addEventListener("click", () => setViewMode(element.dataset.viewMode));
  const requested = new URLSearchParams(location.search).get("run");
  await selectRun(state.runs.some((run) => run.name === requested) ? requested : state.runs[0].name);
}

requestAnimationFrame(animate);
init().catch(showFatal);
