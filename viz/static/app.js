import * as THREE from "/vendor/three.module.js";
import { buildStreamlines, wrapLongitude } from "/streamlines.js";

const $ = (selector) => document.querySelector(selector);
const state = {
  runs: [], metadata: null, run: null, field: null, frameIndex: 0, fieldValues: null,
  frameRecord: null, scaleMode: "auto", scaleMin: 0, scaleMax: 1, lockedScale: null,
  level: null,
  playing: false, rate: 1, lastTick: 0, accumulator: 0,
  frameCache: new Map(), conservation: null, metricIndex: 0,
  globalScaleCache: new Map(), scaleRequestId: 0,
  cacheGeneration: 0, frameRequestId: 0, frameAbortController: null, playbackLoadPending: false,
  streamlinesEnabled: true, streamlineDensity: "medium", windCache: new Map(), windRecord: null,
  windRequestId: 0, windAbortController: null,
  viewMode: "globe",
};

const ui = {
  run: $("#run-select"), field: $("#field-select"), datasetState: $("#dataset-state"), elapsed: $("#elapsed"),
  points: $("#point-count"), maxValue: $("#max-value"), maxValueLabel: $("#max-value-label"), globe: $("#globe"), map: $("#map"),
  play: $("#play"), timeline: $("#timeline"), frame: $("#frame-number"), time: $("#time-label"),
  rate: $("#rate"), legendMin: $("#legend-min"), legendMax: $("#legend-max"), coordinate: $("#globe-coordinate"),
  scaleMode: $("#scale-mode"), gridPoints: $("#grid-points"),
  streamlines: $("#streamlines"), streamlineDensity: $("#streamline-density"),
  flowState: $("#flow-state"), flowKey: $("#flow-key"),
  levelPicker: $("#level-picker"), level: $("#level-select"),
  tabs: $("#metric-tabs"), chart: $("#chart"), metricCurrent: $("#metric-current"),
  metricDrift: $("#metric-drift"), diagnosticsState: $("#diagnostics-state"),
  diagnostics: $("#diagnostics"), fatal: $("#fatal"),
  globePanel: $("#globe-panel"), mapPanel: $("#map-panel"),
  viewButtons: [...document.querySelectorAll("[data-view-mode]")],
};

const views = { globe: null, map: null };

async function fetchJSON(url) {
  const response = await fetch(url, { cache: "no-store" });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `${response.status} ${response.statusText}`);
  return data;
}

function showFatal(error) {
  console.error(error);
  ui.fatal.textContent = `Visualizer error: ${error.message}`;
  ui.fatal.hidden = false;
}

function stopPlayback() {
  state.playing = false;
  state.accumulator = 0;
  ui.play.classList.remove("playing");
  ui.play.setAttribute("aria-label", "再生");
}

function formatDuration(seconds) {
  const total = Math.max(0, Math.round(seconds));
  const days = Math.floor(total / 86400);
  const hours = Math.floor((total % 86400) / 3600);
  const minutes = Math.floor((total % 3600) / 60);
  const secs = total % 60;
  return `${days ? `${days}D ` : ""}${String(hours).padStart(2, "0")}:${String(minutes).padStart(2, "0")}:${String(secs).padStart(2, "0")}`;
}

function formatValue(value) {
  if (!Number.isFinite(value)) return "—";
  const magnitude = Math.abs(value);
  if ((magnitude > 0 && magnitude < 0.001) || magnitude >= 100000) return value.toExponential(4);
  return value.toLocaleString("ja-JP", { maximumFractionDigits: 5 });
}

function runMaximumSpeed(metadata) {
  const maximumCfl = Number(metadata.numerics?.maximum_cfl);
  const truncation = Number(metadata.numerics?.spectral_truncation);
  const timeStep = Number(metadata.simulation?.time_step_seconds);
  const radius = Number(metadata.physical_constants?.earth_radius_m);
  const spectralScale = Math.sqrt(truncation * (truncation + 1));
  const maximumSpeed = maximumCfl * radius / (timeStep * spectralScale);
  return Number.isFinite(maximumSpeed) && maximumSpeed > 0 ? maximumSpeed : null;
}

function selectedField() {
  return state.metadata?.available_fields.find((field) => field.id === state.field) ?? null;
}

function selectedFieldUsesLevel() {
  return Boolean(selectedField()?.uses_level);
}

function formatReferencePressure(pressurePa) {
  const pressureHpa = pressurePa / 100;
  return pressureHpa >= 100
    ? pressureHpa.toFixed(0)
    : pressureHpa.toFixed(1);
}

function formatFieldValue(value, unit) {
  return `${formatValue(value)} ${unit}`;
}

function percentile(values, quantile, absolute = false) {
  const ordered = Array.from(values, (value) => absolute ? Math.abs(value) : value).sort((a, b) => a - b);
  const index = Math.min(ordered.length - 1, Math.max(0, Math.ceil(quantile * ordered.length) - 1));
  return ordered[index];
}

function createRenderer(container) {
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "high-performance" });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 0);
  container.appendChild(renderer.domElement);
  return renderer;
}

function lineMaterial(opacity = 0.18) {
  return new THREE.LineBasicMaterial({ color: 0x9cc6d5, transparent: true, opacity });
}

function addGlobeGraticule(scene) {
  const material = lineMaterial(0.13);
  for (let latitude = -60; latitude <= 60; latitude += 30) {
    const phi = THREE.MathUtils.degToRad(latitude);
    const points = [];
    for (let i = 0; i <= 128; i += 1) {
      const lon = (i / 128) * Math.PI * 2;
      points.push(new THREE.Vector3(Math.cos(phi) * Math.cos(lon), Math.sin(phi), -Math.cos(phi) * Math.sin(lon)).multiplyScalar(1.012));
    }
    scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
  }
  for (let longitude = 0; longitude < 180; longitude += 30) {
    const points = [];
    for (let i = 0; i <= 128; i += 1) {
      const phi = -Math.PI / 2 + (i / 128) * Math.PI;
      const lon = THREE.MathUtils.degToRad(longitude);
      points.push(new THREE.Vector3(Math.cos(phi) * Math.cos(lon), Math.sin(phi), -Math.cos(phi) * Math.sin(lon)).multiplyScalar(1.012));
    }
    const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material);
    scene.add(line, line.clone().rotateY(Math.PI));
  }
}

function addMapGraticule(scene) {
  const vertices = [];
  for (let longitude = -150; longitude <= 150; longitude += 30) {
    const x = longitude / 180;
    vertices.push(x, -1, 0.01, x, 1, 0.01);
  }
  for (let latitude = -60; latitude <= 60; latitude += 30) {
    const y = latitude / 90;
    vertices.push(-1, y, 0.01, 1, y, 0.01);
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(vertices, 3));
  scene.add(new THREE.LineSegments(geometry, lineMaterial(0.14)));
}

function buildCoordinates(metadata) {
  const { mu, nlon } = metadata.grid;
  const count = metadata.grid.point_count;
  const globe = new Float32Array(count * 3);
  const map = new Float32Array(count * 3);
  let point = 0;
  for (let j = 0; j < mu.length; j += 1) {
    const latitude = Math.asin(mu[j]);
    const cosLatitude = Math.cos(latitude);
    for (let k = 0; k < nlon[j]; k += 1) {
      const longitude = (2 * Math.PI * k) / nlon[j] - Math.PI;
      const index = point * 3;
      globe[index] = 1.008 * cosLatitude * Math.cos(longitude);
      globe[index + 1] = 1.008 * Math.sin(latitude);
      // Mirror the WebGL z axis so eastward (increasing longitude) is screen-right,
      // matching the equirectangular map when the globe is viewed from outside.
      globe[index + 2] = -1.008 * cosLatitude * Math.sin(longitude);
      map[index] = longitude / Math.PI;
      map[index + 1] = latitude / (Math.PI / 2);
      map[index + 2] = 0;
      point += 1;
    }
  }
  return { globe, map };
}

function buildRingTriangles(metadata) {
  const { nlon, ring_offsets: offsets } = metadata.grid;
  const triangles = [];
  for (let ring = 0; ring < nlon.length - 1; ring += 1) {
    const lowerCount = nlon[ring];
    const upperCount = nlon[ring + 1];
    const lowerStart = offsets[ring];
    const upperStart = offsets[ring + 1];
    let lower = 0;
    let upper = 0;
    while (lower < lowerCount || upper < upperCount) {
      const advanceLower = upper >= upperCount || (
        lower < lowerCount
        && (lower + 1) * upperCount < (upper + 1) * lowerCount
      );
      if (advanceLower) {
        triangles.push(
          lowerStart + (lower % lowerCount),
          lowerStart + ((lower + 1) % lowerCount),
          upperStart + (upper % upperCount),
        );
        lower += 1;
      } else {
        triangles.push(
          lowerStart + (lower % lowerCount),
          upperStart + ((upper + 1) % upperCount),
          upperStart + (upper % upperCount),
        );
        upper += 1;
      }
    }
  }
  return triangles;
}

function createGlobeSurface(positions, triangles) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(new Float32Array(positions.length), 3));
  geometry.setIndex(triangles);
  const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
  return { mesh: new THREE.Mesh(geometry, material), colorSources: null };
}

function createMapSurface(positions, triangles) {
  const expandedPositions = [];
  const colorSources = [];
  const appendTriangle = (sources, transformX) => {
    for (const source of sources) {
      const offset = source * 3;
      expandedPositions.push(transformX(positions[offset]), positions[offset + 1], 0);
      colorSources.push(source);
    }
  };
  for (let index = 0; index < triangles.length; index += 3) {
    const sources = triangles.slice(index, index + 3);
    const xs = sources.map((source) => positions[source * 3]);
    if (Math.max(...xs) - Math.min(...xs) > 1) {
      appendTriangle(sources, (x) => (x < 0 ? x + 2 : x));
    } else {
      appendTriangle(sources, (x) => x);
    }
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.Float32BufferAttribute(expandedPositions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(new Float32Array(expandedPositions.length), 3));
  const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
  return { mesh: new THREE.Mesh(geometry, material), colorSources: new Uint32Array(colorSources) };
}

function createPointOverlay(positions, isMap) {
  const overlayPositions = positions.slice();
  if (isMap) {
    for (let index = 2; index < overlayPositions.length; index += 3) overlayPositions[index] = 0.02;
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(overlayPositions, 3));
  const material = new THREE.PointsMaterial({
    size: isMap ? 1.8 : 0.012, color: 0xffffff, sizeAttenuation: !isMap,
    transparent: true, opacity: 0.42, depthTest: false, depthWrite: false,
  });
  const points = new THREE.Points(geometry, material);
  points.visible = ui.gridPoints.checked;
  points.renderOrder = 3;
  return points;
}

function clearStreamlineLayers() {
  for (const view of Object.values(views)) {
    if (!view?.streamlineGroup) continue;
    for (const object of [...view.streamlineGroup.children]) {
      view.streamlineGroup.remove(object);
      object.geometry?.dispose();
      object.material?.dispose();
    }
    view.tracer = null;
  }
}

function flowColor(speed, maximumSpeed, target) {
  const amount = Math.sqrt(THREE.MathUtils.clamp(speed / Math.max(maximumSpeed, Number.EPSILON), 0, 1));
  target.push(
    THREE.MathUtils.lerp(0.20, 0.78, amount),
    THREE.MathUtils.lerp(0.60, 0.95, amount),
    THREE.MathUtils.lerp(0.92, 0.44, amount),
  );
}

function appendFlowPosition(target, point, isMap) {
  if (isMap) {
    target.push(wrapLongitude(point.longitude) / Math.PI, point.latitude / (Math.PI / 2), 0.055);
    return;
  }
  const radius = 1.022;
  const cosLatitude = Math.cos(point.latitude);
  target.push(
    radius * cosLatitude * Math.cos(point.longitude),
    radius * Math.sin(point.latitude),
    -radius * cosLatitude * Math.sin(point.longitude),
  );
}

function renderStreamlines(result) {
  clearStreamlineLayers();
  if (!state.streamlinesEnabled || !result?.lines.length) {
    ui.flowKey.hidden = true;
    ui.flowState.textContent = state.streamlinesEnabled ? "STREAMLINES · NO FLOW" : "STREAMLINES · OFF";
    return;
  }
  ui.flowKey.hidden = false;

  for (const [key, view] of Object.entries(views)) {
    const isMap = key === "map";
    const positions = [];
    const colors = [];
    for (const line of result.lines) {
      for (let index = 1; index < line.length; index += 1) {
        const previous = line[index - 1];
        const current = line[index];
        if (isMap && Math.abs(wrapLongitude(current.longitude) - wrapLongitude(previous.longitude)) > Math.PI) continue;
        appendFlowPosition(positions, previous, isMap);
        appendFlowPosition(positions, current, isMap);
        flowColor(previous.speed, result.maximumSpeed, colors);
        flowColor(current.speed, result.maximumSpeed, colors);
      }
    }
    const lineGeometry = new THREE.BufferGeometry();
    lineGeometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    lineGeometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
    const lineObject = new THREE.LineSegments(lineGeometry, new THREE.LineBasicMaterial({
      vertexColors: true, transparent: true, opacity: 0.82,
      depthTest: !isMap, depthWrite: false,
    }));
    lineObject.renderOrder = 6;
    view.streamlineGroup.add(lineObject);

    const tracerPositions = new Float32Array(result.lines.length * 3);
    const tracerGeometry = new THREE.BufferGeometry();
    tracerGeometry.setAttribute("position", new THREE.BufferAttribute(tracerPositions, 3));
    const tracerObject = new THREE.Points(tracerGeometry, new THREE.PointsMaterial({
      color: 0xe8ff9b, size: isMap ? 3.2 : 0.024, sizeAttenuation: !isMap,
      transparent: true, opacity: 0.96, depthTest: !isMap, depthWrite: false,
    }));
    tracerObject.renderOrder = 7;
    view.streamlineGroup.add(tracerObject);
    const speedFactors = result.lines.map((line) => {
      const meanSpeed = line.reduce((sum, point) => sum + point.speed, 0) / line.length;
      return 0.25 + 0.85 * meanSpeed / Math.max(result.maximumSpeed, Number.EPSILON);
    });
    view.tracer = { object: tracerObject, lines: result.lines, speedFactors, isMap };
  }
  const levelLabel = state.metadata?.vertical_coordinate && Number.isInteger(state.level)
    ? ` · L${String(state.level).padStart(2, "0")}`
    : "";
  ui.flowState.textContent = `STREAMLINES${levelLabel} · ${result.lines.length} PATHS`;
}

function updateFlowTracers(timestamp) {
  for (const view of Object.values(views)) {
    const tracer = view?.tracer;
    if (!tracer) continue;
    const positions = tracer.object.geometry.getAttribute("position");
    tracer.lines.forEach((line, lineIndex) => {
      const phase = (timestamp * 0.00004 * tracer.speedFactors[lineIndex] + lineIndex * 0.61803398875) % 1;
      const point = line[Math.min(line.length - 1, Math.floor(phase * line.length))];
      const offset = lineIndex * 3;
      if (tracer.isMap) {
        positions.array[offset] = wrapLongitude(point.longitude) / Math.PI;
        positions.array[offset + 1] = point.latitude / (Math.PI / 2);
        positions.array[offset + 2] = 0.055;
      } else {
        const radius = 1.022;
        const cosLatitude = Math.cos(point.latitude);
        positions.array[offset] = radius * cosLatitude * Math.cos(point.longitude);
        positions.array[offset + 1] = radius * Math.sin(point.latitude);
        positions.array[offset + 2] = -radius * cosLatitude * Math.sin(point.longitude);
      }
    });
    positions.needsUpdate = true;
  }
}

function setViewMode(mode) {
  if (!Object.hasOwn(views, mode)) throw new Error(`Unknown view mode: ${mode}`);
  state.viewMode = mode;
  ui.globePanel.hidden = mode !== "globe";
  ui.mapPanel.hidden = mode !== "map";
  for (const button of ui.viewButtons) {
    const selected = button.dataset.viewMode === mode;
    button.classList.toggle("active", selected);
    button.setAttribute("aria-pressed", String(selected));
  }
  requestAnimationFrame(() => views[mode]?.resize?.());
}

function setupViews(metadata) {
  for (const key of ["globe", "map"]) {
    if (views[key]) {
      views[key].observer.disconnect();
      views[key].renderer.dispose();
      ui[key].replaceChildren();
    }
  }
  const coordinates = buildCoordinates(metadata);
  const triangles = buildRingTriangles(metadata);

  const globeRenderer = createRenderer(ui.globe);
  const globeScene = new THREE.Scene();
  const globeCamera = new THREE.PerspectiveCamera(38, 1, 0.1, 20);
  globeCamera.position.set(0, 0, 3.25);
  const sphere = new THREE.Mesh(
    new THREE.SphereGeometry(1, 64, 40),
    new THREE.MeshBasicMaterial({ color: 0x071b29, depthWrite: true }),
  );
  globeScene.add(sphere);
  addGlobeGraticule(globeScene);
  const globeSurface = createGlobeSurface(coordinates.globe, triangles);
  const globePoints = createPointOverlay(coordinates.globe, false);
  const globeStreamlines = new THREE.Group();
  globeScene.add(globeSurface.mesh, globePoints, globeStreamlines);
  views.globe = {
    renderer: globeRenderer, scene: globeScene, camera: globeCamera,
    surface: globeSurface.mesh, colorSources: globeSurface.colorSources,
    pointOverlay: globePoints, streamlineGroup: globeStreamlines, tracer: null,
    yaw: Math.PI / 2, pitch: 0.28, distance: 3.25,
  };
  installGlobeControls(views.globe);

  const mapRenderer = createRenderer(ui.map);
  const mapScene = new THREE.Scene();
  const mapCamera = new THREE.OrthographicCamera(-1.08, 1.08, 1.08, -1.08, 0.1, 10);
  mapCamera.position.z = 2;
  addMapGraticule(mapScene);
  const mapSurface = createMapSurface(coordinates.map, triangles);
  const mapPoints = createPointOverlay(coordinates.map, true);
  const mapStreamlines = new THREE.Group();
  mapScene.add(mapSurface.mesh, mapPoints, mapStreamlines);
  views.map = {
    renderer: mapRenderer, scene: mapScene, camera: mapCamera,
    surface: mapSurface.mesh, colorSources: mapSurface.colorSources,
    pointOverlay: mapPoints, streamlineGroup: mapStreamlines, tracer: null,
  };

  for (const key of ["globe", "map"]) {
    const view = views[key];
    const resize = () => {
      const width = ui[key].clientWidth;
      const height = ui[key].clientHeight;
      view.renderer.setSize(Math.max(1, width), Math.max(1, height), false);
      if (key === "globe") {
        view.camera.aspect = width / Math.max(1, height);
      } else {
        // Both axes fill the viewport; the CSS viewport itself is exactly 2:1,
        // so 360 degrees of longitude occupy twice the pixels of 180 degrees latitude.
        view.camera.left = -1.08;
        view.camera.right = 1.08;
        view.camera.top = 1.08;
        view.camera.bottom = -1.08;
      }
      view.camera.updateProjectionMatrix();
    };
    view.resize = resize;
    view.observer = new ResizeObserver(resize);
    view.observer.observe(ui[key]);
    resize();
  }
  setViewMode(state.viewMode);
}

function installGlobeControls(view) {
  let dragging = false;
  let previous = { x: 0, y: 0 };
  const canvas = view.renderer.domElement;
  canvas.addEventListener("pointerdown", (event) => {
    dragging = true; previous = { x: event.clientX, y: event.clientY }; canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener("pointermove", (event) => {
    if (dragging) {
      view.yaw += (event.clientX - previous.x) * 0.006;
      view.pitch = THREE.MathUtils.clamp(view.pitch + (event.clientY - previous.y) * 0.006, -1.35, 1.35);
      previous = { x: event.clientX, y: event.clientY };
    }
    const rect = canvas.getBoundingClientRect();
    const lon = ((event.clientX - rect.left) / rect.width) * 360 - 180;
    const lat = 90 - ((event.clientY - rect.top) / rect.height) * 180;
    ui.coordinate.textContent = `λ ${lon >= 0 ? "+" : "−"}${Math.abs(lon).toFixed(1)}° · φ ${lat >= 0 ? "+" : "−"}${Math.abs(lat).toFixed(1)}°`;
  });
  canvas.addEventListener("pointerup", () => { dragging = false; });
  canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    view.distance = THREE.MathUtils.clamp(view.distance + event.deltaY * 0.002, 2.25, 5.2);
  }, { passive: false });
}

function colorFor(t, target, offset) {
  const stops = [
    [0.00, [0.082, 0.369, 0.937]], [0.28, [0.192, 0.710, 0.925]],
    [0.53, [0.898, 0.933, 0.949]], [0.76, [1.000, 0.624, 0.263]], [1.00, [0.941, 0.267, 0.220]],
  ];
  let upper = 1;
  while (upper < stops.length - 1 && t > stops[upper][0]) upper += 1;
  const lower = upper - 1;
  const amount = (t - stops[lower][0]) / (stops[upper][0] - stops[lower][0]);
  for (let channel = 0; channel < 3; channel += 1) {
    target[offset + channel] = THREE.MathUtils.lerp(stops[lower][1][channel], stops[upper][1][channel], amount);
  }
}

function updateColors(values, scaleMin, scaleMax) {
  const scaleWidth = Math.max(scaleMax - scaleMin, Number.EPSILON);
  const pointColors = new Float32Array(values.length * 3);
  for (let i = 0; i < values.length; i += 1) {
    if (!Number.isFinite(values[i])) {
      pointColors.set([0.12, 0.14, 0.17], i * 3);
      continue;
    }
    const normalized = (values[i] - scaleMin) / scaleWidth;
    colorFor(THREE.MathUtils.clamp(normalized, 0, 1), pointColors, i * 3);
  }
  for (const key of ["globe", "map"]) {
    const view = views[key];
    const attribute = view.surface.geometry.getAttribute("color");
    if (view.colorSources === null) {
      attribute.array.set(pointColors);
    } else {
      for (let vertex = 0; vertex < view.colorSources.length; vertex += 1) {
        const sourceOffset = view.colorSources[vertex] * 3;
        const targetOffset = vertex * 3;
        attribute.array[targetOffset] = pointColors[sourceOffset];
        attribute.array[targetOffset + 1] = pointColors[sourceOffset + 1];
        attribute.array[targetOffset + 2] = pointColors[sourceOffset + 2];
      }
    }
    attribute.needsUpdate = true;
  }
}

function renderFieldRecord(record) {
  const field = selectedField();
  if (!field || !record) return;
  if (state.scaleMode === "auto") {
    if (field.signed) {
      state.scaleMax = Math.max(record.p995Absolute, Number.EPSILON);
      state.scaleMin = -state.scaleMax;
    } else if (field.zero_based) {
      state.scaleMin = 0;
      state.scaleMax = Math.max(record.p995, Number.EPSILON);
    } else {
      state.scaleMin = record.p005;
      state.scaleMax = Math.max(record.p995, state.scaleMin + Number.EPSILON);
    }
  } else if (state.scaleMode === "lock") {
    state.scaleMin = state.lockedScale?.minimum ?? state.scaleMin;
    state.scaleMax = state.lockedScale?.maximum ?? state.scaleMax;
  } else if (!(state.scaleMax > state.scaleMin)) {
    state.scaleMin = field.signed ? -Math.max(record.p995Absolute, Number.EPSILON) : 0;
    state.scaleMax = Math.max(record.p995Absolute, Number.EPSILON);
  }
  updateColors(record.values, state.scaleMin, state.scaleMax);
  const frameMaximum = field.signed
    ? Math.max(Math.abs(record.minimum), Math.abs(record.maximum))
    : record.maximum;
  ui.maxValue.textContent = formatFieldValue(frameMaximum, field.unit);
  ui.legendMin.textContent = formatFieldValue(state.scaleMin, field.unit);
  ui.legendMax.textContent = formatFieldValue(state.scaleMax, field.unit);
}

async function setScaleMode(mode) {
  if (!["auto", "global", "lock"].includes(mode)) throw new Error(`Unknown scale mode: ${mode}`);
  state.scaleMode = mode;
  ui.scaleMode.value = mode;
  const requestId = ++state.scaleRequestId;
  if (mode === "lock") {
    state.lockedScale = { minimum: state.scaleMin, maximum: state.scaleMax };
    renderFieldRecord(state.frameRecord);
    return;
  }
  state.lockedScale = null;
  if (mode === "auto") {
    renderFieldRecord(state.frameRecord);
    return;
  }

  const runAtStart = state.run;
  const fieldAtStart = state.field;
  const levelAtStart = state.level;
  const levelKey = selectedFieldUsesLevel() ? state.level : "surface";
  const cacheKey = `${runAtStart}:${fieldAtStart}:${levelKey}`;
  let scale = state.globalScaleCache.get(cacheKey);
  if (!scale && fieldAtStart === "speed") {
    const speedMaximum = runMaximumSpeed(state.metadata);
    if (speedMaximum > 0) scale = { minimum: 0, maximum: speedMaximum };
  }
  if (!scale) {
    const query = selectedFieldUsesLevel() && Number.isInteger(levelAtStart)
      ? `?level=${encodeURIComponent(levelAtStart)}`
      : "";
    const stats = await fetchJSON(`/api/runs/${encodeURIComponent(runAtStart)}/fields/${encodeURIComponent(fieldAtStart)}/statistics${query}`);
    const field = selectedField();
    scale = field?.signed
      ? { minimum: -stats.maximum_absolute, maximum: stats.maximum_absolute }
      : { minimum: field?.zero_based ? 0 : stats.minimum, maximum: stats.maximum };
  }
  if (
    requestId !== state.scaleRequestId
    || state.run !== runAtStart
    || state.field !== fieldAtStart
    || state.level !== levelAtStart
    || state.scaleMode !== "global"
  ) return;
  state.scaleMin = scale.minimum;
  state.scaleMax = Math.max(scale.maximum, scale.minimum + Number.EPSILON);
  state.globalScaleCache.set(cacheKey, { minimum: state.scaleMin, maximum: state.scaleMax });
  renderFieldRecord(state.frameRecord);
}

function renderWindRecord(record) {
  if (!record || !state.metadata) return;
  renderStreamlines(buildStreamlines(
    state.metadata.grid,
    record.eastward,
    record.northward,
    state.streamlineDensity,
  ));
}

async function loadWindFrame(frameIndex) {
  if (!state.streamlinesEnabled || !state.metadata?.supports_streamlines) {
    clearStreamlineLayers();
    ui.flowKey.hidden = true;
    ui.flowState.textContent = state.metadata?.supports_streamlines === false
      ? "STREAMLINES · UNAVAILABLE"
      : "STREAMLINES · OFF";
    return false;
  }
  const runAtStart = state.run;
  const metadataAtStart = state.metadata;
  const levelAtStart = state.level;
  const generationAtStart = state.cacheGeneration;
  const step = metadataAtStart.available_steps[frameIndex];
  const levelKey = metadataAtStart.vertical_coordinate ? levelAtStart : "surface";
  const cacheKey = `${generationAtStart}:${runAtStart}:${levelKey}:${step}`;
  const requestId = ++state.windRequestId;
  if (state.windAbortController) state.windAbortController.abort();
  state.windAbortController = null;
  ui.flowState.textContent = "STREAMLINES · TRACING…";

  let record = state.windCache.get(cacheKey);
  if (!record) {
    const controller = new AbortController();
    state.windAbortController = controller;
    try {
      const params = new URLSearchParams({ generation: String(generationAtStart) });
      if (metadataAtStart.vertical_coordinate && Number.isInteger(levelAtStart)) {
        params.set("level", String(levelAtStart));
      }
      const response = await fetch(
        `/api/runs/${encodeURIComponent(runAtStart)}/wind/${step}?${params}`,
        { signal: controller.signal, cache: "no-store" },
      );
      if (!response.ok) {
        const error = await response.json();
        throw new Error(error.error || "Could not load wind components");
      }
      const values = new Float32Array(await response.arrayBuffer());
      const count = metadataAtStart.grid.point_count;
      if (values.length !== count * 2) throw new Error("Wind component point count does not match metadata");
      record = {
        eastward: values.subarray(0, count),
        northward: values.subarray(count),
      };
      state.windCache.set(cacheKey, record);
      while (state.windCache.size > 8) state.windCache.delete(state.windCache.keys().next().value);
    } catch (error) {
      if (error.name === "AbortError") return false;
      if (requestId === state.windRequestId) {
        clearStreamlineLayers();
        ui.flowKey.hidden = true;
        ui.flowState.textContent = "STREAMLINES · ERROR";
      }
      console.error(error);
      return false;
    } finally {
      if (state.windAbortController === controller) state.windAbortController = null;
    }
  }
  if (
    requestId !== state.windRequestId
    || state.run !== runAtStart
    || state.metadata !== metadataAtStart
    || state.level !== levelAtStart
    || state.frameIndex !== frameIndex
    || state.cacheGeneration !== generationAtStart
    || !state.streamlinesEnabled
  ) return false;
  state.windRecord = record;
  renderWindRecord(record);
  return true;
}

async function loadFrame(frameIndex) {
  const runAtStart = state.run;
  const fieldAtStart = state.field;
  const metadataAtStart = state.metadata;
  const levelAtStart = state.level;
  const cacheGenerationAtStart = state.cacheGeneration;
  const requestId = ++state.frameRequestId;
  if (state.frameAbortController) state.frameAbortController.abort();
  state.frameAbortController = null;

  const step = metadataAtStart.available_steps[frameIndex];
  const fieldAtStartMetadata = metadataAtStart.available_fields.find((field) => field.id === fieldAtStart);
  const levelKey = fieldAtStartMetadata?.uses_level ? levelAtStart : "surface";
  const key = `${cacheGenerationAtStart}:${runAtStart}:${fieldAtStart}:${levelKey}:${step}`;
  let record = state.frameCache.get(key);
  if (!record) {
    const controller = new AbortController();
    state.frameAbortController = controller;
    try {
      const params = new URLSearchParams({ generation: String(cacheGenerationAtStart) });
      if (fieldAtStartMetadata?.uses_level && Number.isInteger(levelAtStart)) {
        params.set("level", String(levelAtStart));
      }
      const response = await fetch(
        `/api/runs/${encodeURIComponent(runAtStart)}/fields/${encodeURIComponent(fieldAtStart)}/${step}?${params}`,
        { signal: controller.signal, cache: "no-store" },
      );
      if (!response.ok) {
        const error = await response.json();
        if (requestId !== state.frameRequestId || state.run !== runAtStart) return false;
        throw new Error(error.error || "Could not load frame");
      }
      const values = new Float32Array(await response.arrayBuffer());
      if (values.length !== metadataAtStart.grid.point_count) throw new Error("Frame point count does not match metadata");
      const percentileHeader = response.headers.get("X-Field-P995-Absolute");
      const headerPercentile = percentileHeader === null ? Number.NaN : Number(percentileHeader);
      const p005Header = response.headers.get("X-Field-P005");
      const p995Header = response.headers.get("X-Field-P995");
      record = {
        values,
        minimum: Number(response.headers.get("X-Field-Minimum")),
        maximum: Number(response.headers.get("X-Field-Maximum")),
        p005: p005Header === null ? Number.NaN : Number(p005Header),
        p995: p995Header === null ? Number.NaN : Number(p995Header),
        p995Absolute: Number.isFinite(headerPercentile) && headerPercentile >= 0
          ? headerPercentile
          : percentile(values, 0.995, true),
      };
      if (!Number.isFinite(record.p005)) record.p005 = percentile(values, 0.005);
      if (!Number.isFinite(record.p995)) record.p995 = percentile(values, 0.995);
      state.frameCache.set(key, record);
      while (state.frameCache.size > 16) state.frameCache.delete(state.frameCache.keys().next().value);
    } catch (error) {
      if (error.name === "AbortError") return false;
      throw error;
    } finally {
      if (state.frameAbortController === controller) state.frameAbortController = null;
    }
  }
  if (
    state.run !== runAtStart
    || state.field !== fieldAtStart
    || state.metadata !== metadataAtStart
    || state.level !== levelAtStart
    || state.cacheGeneration !== cacheGenerationAtStart
    || requestId !== state.frameRequestId
  ) return false;
  state.frameIndex = frameIndex;
  state.fieldValues = record.values;
  state.frameRecord = record;
  renderFieldRecord(record);
  ui.timeline.value = String(frameIndex);
  ui.frame.textContent = String(step).padStart(5, "0");
  const seconds = state.metadata.frame_times_seconds?.[String(step)]
    ?? step * state.metadata.simulation.time_step_seconds;
  const sampling = state.metadata.surface_sampling;
  ui.time.textContent = `${sampling === "monthly" ? "MONTH MEAN · " : sampling === "yearly" ? "SNAPSHOT · " : ""}T + ${formatDuration(seconds)}`;
  ui.elapsed.textContent = formatDuration(seconds);
  updateChartReadout();
  drawChart();
  await loadWindFrame(frameIndex);
  return true;
}

async function setField(fieldId, frameIndex = state.frameIndex) {
  stopPlayback();
  state.frameRequestId += 1;
  if (state.frameAbortController) state.frameAbortController.abort();
  state.frameAbortController = null;
  state.playbackLoadPending = false;
  state.scaleRequestId += 1;
  const field = state.metadata.available_fields.find((candidate) => candidate.id === fieldId);
  if (!field) throw new Error(`Field ${fieldId} is not available for this run`);
  state.field = fieldId;
  state.frameRecord = null;
  state.scaleMin = 0;
  state.scaleMax = 0;
  state.lockedScale = null;
  if (state.scaleMode === "lock") {
    state.scaleMode = "auto";
    ui.scaleMode.value = "auto";
  }
  ui.field.value = fieldId;
  ui.maxValueLabel.textContent = field.signed ? `MAX |${field.symbol}|` : `MAX ${field.symbol}`;
  ui.globe.setAttribute("aria-label", `${field.label} の三次元地球表示`);
  ui.map.setAttribute("aria-label", `${field.label} の二次元正距円筒図法表示`);
  ui.level.disabled = !field.uses_level && !state.streamlinesEnabled;

  await loadFrame(frameIndex);
  if (state.field === fieldId && state.scaleMode === "global") await setScaleMode("global");
}

async function setLevelFromInput() {
  if (!state.metadata?.vertical_coordinate) return;
  const vertical = state.metadata.vertical_coordinate;
  const requested = Number(ui.level.value);
  if (!Number.isInteger(requested)) {
    ui.level.value = String(state.level);
    return;
  }
  const level = Math.min(
    Number(vertical.number_of_levels),
    Math.max(1, Math.trunc(requested)),
  );
  ui.level.value = String(level);
  if (level === state.level) return;
  state.level = level;
  stopPlayback();
  if (!selectedFieldUsesLevel()) {
    await loadWindFrame(state.frameIndex);
    return;
  }
  state.frameRecord = null;
  state.scaleMin = 0;
  state.scaleMax = 0;
  state.lockedScale = null;
  if (state.scaleMode === "lock") {
    state.scaleMode = "auto";
    ui.scaleMode.value = "auto";
  }
  await loadFrame(state.frameIndex);
  if (state.scaleMode === "global") await setScaleMode("global");
}

async function selectRun(name) {
  stopPlayback();
  state.cacheGeneration += 1;
  state.frameCache.clear();
  state.windCache.clear();
  state.windRecord = null;
  state.globalScaleCache.clear();
  state.frameRequestId += 1;
  if (state.frameAbortController) state.frameAbortController.abort();
  state.frameAbortController = null;
  state.windRequestId += 1;
  if (state.windAbortController) state.windAbortController.abort();
  state.windAbortController = null;
  state.playbackLoadPending = false;
  state.run = name;
  state.field = null;
  state.level = null;
  state.frameRecord = null;
  state.scaleRequestId += 1;
  state.frameIndex = 0;
  state.conservation = null;
  state.metricIndex = 0;
  ui.diagnosticsState.textContent = "COMPUTING SPHERICAL INTEGRALS…";
  ui.tabs.replaceChildren();
  const metadata = await fetchJSON(`/api/runs/${encodeURIComponent(name)}/metadata`);
  if (state.run !== name) return;
  state.metadata = metadata;
  ui.streamlines.disabled = !metadata.supports_streamlines;
  ui.streamlineDensity.disabled = !metadata.supports_streamlines || !state.streamlinesEnabled;
  const vertical = metadata.vertical_coordinate;
  ui.levelPicker.hidden = !vertical;
  if (vertical) {
    state.level = Number(vertical.default_level);
    ui.level.replaceChildren(...vertical.reference_full_level_pressure_pa.map((pressurePa, index) => {
      const option = document.createElement("option");
      option.value = String(index + 1);
      option.textContent = `L${String(index + 1).padStart(2, "0")} · ≈ ${formatReferencePressure(Number(pressurePa))} hPa`;
      return option;
    }));
    ui.level.value = String(state.level);
  }
  ui.field.replaceChildren(...metadata.available_fields.map((field) => {
    const option = document.createElement("option");
    option.value = field.id;
    option.textContent = field.label;
    return option;
  }));
  ui.points.textContent = Number(metadata.grid.point_count).toLocaleString("ja-JP");
  ui.timeline.max = String(metadata.available_steps.length - 1);
  ui.datasetState.textContent = vertical
    ? `${metadata.available_frame_count} FRAMES · LEVELS READY`
    : `${metadata.available_frame_count} FRAMES READY`;
  setupViews(metadata);
  const initialField = metadata.available_fields.some((field) => field.id === "speed")
    ? "speed"
    : metadata.available_fields[0].id;
  await setField(initialField, 0);
  if (state.run !== name) return;
  ui.diagnostics.hidden = metadata.supports_conservation_diagnostics === false;
  if (metadata.supports_conservation_diagnostics === false) return;
  fetchJSON(`/api/runs/${encodeURIComponent(name)}/conservation`)
    .then((data) => {
      if (state.run !== name) return;
      state.conservation = data;
      ui.diagnosticsState.textContent = "SPHERICAL GAUSSIAN QUADRATURE";
      buildMetricTabs();
      drawChart();
    })
    .catch(showFatal);
}

async function advancePlaybackFrame() {
  if (state.playbackLoadPending || !state.playing || !state.metadata) return;
  state.playbackLoadPending = true;
  const next = (state.frameIndex + 1) % state.metadata.available_steps.length;
  try {
    await loadFrame(next);
  } catch (error) {
    stopPlayback();
    showFatal(error);
  } finally {
    state.playbackLoadPending = false;
  }
}

function buildMetricTabs() {
  ui.tabs.replaceChildren();
  state.conservation.metrics.forEach((metric, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = metric.label;
    button.classList.toggle("active", index === state.metricIndex);
    button.addEventListener("click", () => {
      state.metricIndex = index;
      [...ui.tabs.children].forEach((node, i) => node.classList.toggle("active", i === index));
      updateChartReadout();
      drawChart();
    });
    ui.tabs.appendChild(button);
  });
  updateChartReadout();
}

function updateChartReadout() {
  if (!state.conservation) { ui.metricCurrent.textContent = "—"; ui.metricDrift.textContent = "—"; return; }
  const metric = state.conservation.metrics[state.metricIndex];
  const index = Math.min(state.frameIndex, metric.values.length - 1);
  const current = metric.values[index];
  const initial = metric.values[0];
  ui.metricCurrent.textContent = `${formatValue(current)} ${metric.unit}`;
  if (metric.id === "mean_relative_vorticity") {
    const delta = current - initial;
    ui.metricDrift.textContent = `Δ₀ ${delta >= 0 ? "+" : ""}${delta.toExponential(2)} ${metric.unit}`;
  } else {
    const scale = Math.max(Math.abs(initial), Number.EPSILON);
    const drift = ((current - initial) / scale) * 100;
    ui.metricDrift.textContent = `Δ₀ ${drift >= 0 ? "+" : ""}${drift.toExponential(2)} %`;
  }
}

function drawChart() {
  const canvas = ui.chart;
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.min(devicePixelRatio, 2);
  const width = Math.max(1, Math.round(rect.width * dpr));
  const height = Math.max(1, Math.round(rect.height * dpr));
  if (canvas.width !== width || canvas.height !== height) { canvas.width = width; canvas.height = height; }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, rect.width, rect.height);
  const pad = { left: 62, right: 22, top: 28, bottom: 32 };
  const w = rect.width - pad.left - pad.right;
  const h = rect.height - pad.top - pad.bottom;
  ctx.strokeStyle = "rgba(170, 199, 213, .13)";
  ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i += 1) {
    const y = pad.top + (h * i) / 4;
    ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(pad.left + w, y); ctx.stroke();
  }
  if (!state.conservation) return;
  const metric = state.conservation.metrics[state.metricIndex];
  const values = metric.values;
  let min = Math.min(...values), max = Math.max(...values);
  if (min === max) { min -= Math.abs(min || 1) * 0.01; max += Math.abs(max || 1) * 0.01; }
  const margin = (max - min) * 0.12;
  min -= margin; max += margin;
  ctx.fillStyle = "#6f8492";
  ctx.font = "9px ui-monospace, monospace";
  ctx.textAlign = "right";
  for (let i = 0; i <= 4; i += 1) {
    const value = max - ((max - min) * i) / 4;
    ctx.fillText(Number(value).toExponential(2), pad.left - 9, pad.top + (h * i) / 4 + 3);
  }
  ctx.textAlign = "left";
  ctx.fillText("0 h", pad.left, pad.top + h + 22);
  ctx.textAlign = "right";
  const finalHour = state.conservation.times_seconds.at(-1) / 3600;
  ctx.fillText(`${formatValue(finalHour)} h`, pad.left + w, pad.top + h + 22);
  const gradient = ctx.createLinearGradient(pad.left, 0, pad.left + w, 0);
  gradient.addColorStop(0, "#2388ff"); gradient.addColorStop(1, "#c8f06b");
  ctx.strokeStyle = gradient; ctx.lineWidth = 2;
  ctx.beginPath();
  values.forEach((value, index) => {
    const x = pad.left + (index / Math.max(1, values.length - 1)) * w;
    const y = pad.top + ((max - value) / (max - min)) * h;
    if (index === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  });
  ctx.stroke();
  const markerIndex = Math.min(state.frameIndex, values.length - 1);
  const markerX = pad.left + (markerIndex / Math.max(1, values.length - 1)) * w;
  const markerY = pad.top + ((max - values[markerIndex]) / (max - min)) * h;
  ctx.strokeStyle = "rgba(200,240,107,.4)"; ctx.beginPath(); ctx.moveTo(markerX, pad.top); ctx.lineTo(markerX, pad.top + h); ctx.stroke();
  ctx.fillStyle = "#c8f06b"; ctx.beginPath(); ctx.arc(markerX, markerY, 4, 0, Math.PI * 2); ctx.fill();
}

function animate(timestamp) {
  requestAnimationFrame(animate);
  if (views.globe) {
    const view = views.globe;
    view.camera.position.set(
      view.distance * Math.sin(view.yaw) * Math.cos(view.pitch),
      view.distance * Math.sin(view.pitch),
      view.distance * Math.cos(view.yaw) * Math.cos(view.pitch),
    );
    view.camera.lookAt(0, 0, 0);
    updateFlowTracers(timestamp);
    const activeView = views[state.viewMode];
    activeView.renderer.render(activeView.scene, activeView.camera);
  }
  if (state.playing && state.metadata) {
    const delta = timestamp - state.lastTick;
    state.accumulator += delta * state.rate;
    if (state.accumulator >= 180) {
      state.accumulator %= 180;
      void advancePlaybackFrame();
    }
  }
  state.lastTick = timestamp;
}

async function init() {
  const data = await fetchJSON("/api/runs");
  state.runs = data.runs;
  if (!state.runs.length) throw new Error("output/ に可視化できる計算結果がありません");
  ui.run.replaceChildren(...state.runs.map((run) => {
    const option = document.createElement("option");
    option.value = run.name; option.textContent = run.case_name; return option;
  }));
  ui.run.addEventListener("change", () => selectRun(ui.run.value).catch(showFatal));
  ui.field.addEventListener("change", () => setField(ui.field.value).catch(showFatal));
  ui.level.addEventListener("change", () => setLevelFromInput().catch(showFatal));
  ui.play.addEventListener("click", () => {
    if (state.playing) {
      stopPlayback();
    } else {
      state.playing = true;
      state.accumulator = 0;
      ui.play.classList.add("playing");
      ui.play.setAttribute("aria-label", "一時停止");
    }
  });
  ui.timeline.addEventListener("input", () => {
    stopPlayback();
    loadFrame(Number(ui.timeline.value)).catch(showFatal);
  });
  ui.rate.addEventListener("change", () => { state.rate = Number(ui.rate.value); });
  ui.scaleMode.addEventListener("change", () => setScaleMode(ui.scaleMode.value).catch(showFatal));
  ui.gridPoints.addEventListener("change", () => {
    for (const view of Object.values(views)) {
      if (view) view.pointOverlay.visible = ui.gridPoints.checked;
    }
  });
  ui.streamlines.addEventListener("change", () => {
    state.streamlinesEnabled = ui.streamlines.checked;
    ui.streamlineDensity.disabled = !state.streamlinesEnabled || !state.metadata?.supports_streamlines;
    ui.level.disabled = !selectedFieldUsesLevel() && !state.streamlinesEnabled;
    state.windRequestId += 1;
    if (state.windAbortController) state.windAbortController.abort();
    state.windAbortController = null;
    if (state.streamlinesEnabled) {
      loadWindFrame(state.frameIndex).catch(showFatal);
    } else {
      clearStreamlineLayers();
      ui.flowKey.hidden = true;
      ui.flowState.textContent = "STREAMLINES · OFF";
    }
  });
  ui.streamlineDensity.addEventListener("change", () => {
    state.streamlineDensity = ui.streamlineDensity.value;
    if (state.windRecord && state.streamlinesEnabled) renderWindRecord(state.windRecord);
  });
  for (const button of ui.viewButtons) {
    button.addEventListener("click", () => setViewMode(button.dataset.viewMode));
  }
  window.addEventListener("resize", drawChart);
  await selectRun(state.runs[0].name);
}

requestAnimationFrame(animate);
init().catch(showFatal);
