import * as THREE from "/vendor/three.module.js";

const $ = (selector) => document.querySelector(selector);
const state = {
  runs: [], metadata: null, run: null, frameIndex: 0, speed: null,
  scaleMax: 1, playing: false, rate: 1, lastTick: 0, accumulator: 0,
  frameCache: new Map(), conservation: null, metricIndex: 0,
  frameRequestId: 0, frameAbortController: null, playbackLoadPending: false,
};

const ui = {
  run: $("#run-select"), datasetState: $("#dataset-state"), elapsed: $("#elapsed"),
  points: $("#point-count"), maxSpeed: $("#max-speed"), globe: $("#globe"), map: $("#map"),
  play: $("#play"), timeline: $("#timeline"), frame: $("#frame-number"), time: $("#time-label"),
  rate: $("#rate"), legendMax: $("#legend-max"), coordinate: $("#globe-coordinate"),
  tabs: $("#metric-tabs"), chart: $("#chart"), metricCurrent: $("#metric-current"),
  metricDrift: $("#metric-drift"), diagnosticsState: $("#diagnostics-state"), fatal: $("#fatal"),
};

const views = { globe: null, map: null };

async function fetchJSON(url) {
  const response = await fetch(url);
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
      points.push(new THREE.Vector3(Math.cos(phi) * Math.cos(lon), Math.sin(phi), Math.cos(phi) * Math.sin(lon)).multiplyScalar(1.003));
    }
    scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
  }
  for (let longitude = 0; longitude < 180; longitude += 30) {
    const points = [];
    for (let i = 0; i <= 128; i += 1) {
      const phi = -Math.PI / 2 + (i / 128) * Math.PI;
      const lon = THREE.MathUtils.degToRad(longitude);
      points.push(new THREE.Vector3(Math.cos(phi) * Math.cos(lon), Math.sin(phi), Math.cos(phi) * Math.sin(lon)).multiplyScalar(1.003));
    }
    const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material);
    scene.add(line, line.clone().rotateY(Math.PI));
  }
}

function addMapGraticule(scene) {
  const vertices = [];
  for (let longitude = -150; longitude <= 150; longitude += 30) {
    const x = longitude / 180;
    vertices.push(x, -1, 0, x, 1, 0);
  }
  for (let latitude = -60; latitude <= 60; latitude += 30) {
    const y = latitude / 90;
    vertices.push(-1, y, 0, 1, y, 0);
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
      globe[index] = 1.009 * cosLatitude * Math.cos(longitude);
      globe[index + 1] = 1.009 * Math.sin(latitude);
      globe[index + 2] = 1.009 * cosLatitude * Math.sin(longitude);
      map[index] = longitude / Math.PI;
      map[index + 1] = latitude / (Math.PI / 2);
      map[index + 2] = 0;
      point += 1;
    }
  }
  return { globe, map };
}

function createPointCloud(positions, isMap) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(new Float32Array(positions.length), 3));
  const material = new THREE.PointsMaterial({
    size: isMap ? 2.0 : 0.019, vertexColors: true, sizeAttenuation: !isMap,
    transparent: true, opacity: 0.94, depthWrite: !isMap,
  });
  return new THREE.Points(geometry, material);
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
  const globePoints = createPointCloud(coordinates.globe, false);
  globeScene.add(globePoints);
  views.globe = { renderer: globeRenderer, scene: globeScene, camera: globeCamera, points: globePoints, yaw: -0.5, pitch: 0.28, distance: 3.25 };
  installGlobeControls(views.globe);

  const mapRenderer = createRenderer(ui.map);
  const mapScene = new THREE.Scene();
  const mapCamera = new THREE.OrthographicCamera(-1.08, 1.08, 1.08, -1.08, 0.1, 10);
  mapCamera.position.z = 2;
  addMapGraticule(mapScene);
  const mapPoints = createPointCloud(coordinates.map, true);
  mapScene.add(mapPoints);
  views.map = { renderer: mapRenderer, scene: mapScene, camera: mapCamera, points: mapPoints };

  for (const key of ["globe", "map"]) {
    const view = views[key];
    const resize = () => {
      const width = ui[key].clientWidth;
      const height = ui[key].clientHeight;
      view.renderer.setSize(width, height, false);
      if (key === "globe") {
        view.camera.aspect = width / Math.max(1, height);
      } else {
        const aspect = width / Math.max(1, height);
        view.camera.left = -1.08 * aspect;
        view.camera.right = 1.08 * aspect;
        view.camera.top = 1.08;
        view.camera.bottom = -1.08;
      }
      view.camera.updateProjectionMatrix();
    };
    view.observer = new ResizeObserver(resize);
    view.observer.observe(ui[key]);
    resize();
  }
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

function updateColors(speed, scaleMax) {
  for (const key of ["globe", "map"]) {
    const attribute = views[key].points.geometry.getAttribute("color");
    for (let i = 0; i < speed.length; i += 1) colorFor(Math.min(1, speed[i] / scaleMax), attribute.array, i * 3);
    attribute.needsUpdate = true;
  }
}

async function loadFrame(frameIndex) {
  const runAtStart = state.run;
  const metadataAtStart = state.metadata;
  const requestId = ++state.frameRequestId;
  if (state.frameAbortController) state.frameAbortController.abort();
  state.frameAbortController = null;

  const step = metadataAtStart.available_steps[frameIndex];
  const key = `${runAtStart}:${step}`;
  let record = state.frameCache.get(key);
  if (!record) {
    const controller = new AbortController();
    state.frameAbortController = controller;
    try {
      const response = await fetch(`/api/runs/${encodeURIComponent(runAtStart)}/frame/${step}`, { signal: controller.signal });
      if (!response.ok) {
        const error = await response.json();
        if (requestId !== state.frameRequestId || state.run !== runAtStart) return false;
        throw new Error(error.error || "Could not load frame");
      }
      const speed = new Float32Array(await response.arrayBuffer());
      if (speed.length !== metadataAtStart.grid.point_count) throw new Error("Frame point count does not match metadata");
      record = {
        speed,
        max: Number(response.headers.get("X-Speed-Max")),
        p98: Number(response.headers.get("X-Speed-P98")),
      };
      state.frameCache.set(key, record);
      while (state.frameCache.size > 16) state.frameCache.delete(state.frameCache.keys().next().value);
    } catch (error) {
      if (error.name === "AbortError") return false;
      throw error;
    } finally {
      if (state.frameAbortController === controller) state.frameAbortController = null;
    }
  }
  if (state.run !== runAtStart || state.metadata !== metadataAtStart || requestId !== state.frameRequestId) return false;
  state.frameIndex = frameIndex;
  state.speed = record.speed;
  if (!(state.scaleMax > 0)) state.scaleMax = Math.max(record.p98, Number.EPSILON);
  updateColors(record.speed, state.scaleMax);
  ui.timeline.value = String(frameIndex);
  ui.frame.textContent = String(step).padStart(5, "0");
  const seconds = step * state.metadata.simulation.time_step_seconds;
  ui.time.textContent = `T + ${formatDuration(seconds)}`;
  ui.elapsed.textContent = formatDuration(seconds);
  ui.maxSpeed.textContent = `${record.max.toFixed(2)} m/s`;
  ui.legendMax.textContent = `${state.scaleMax.toFixed(1)} m s⁻¹`;
  updateChartReadout();
  drawChart();
  return true;
}

async function selectRun(name) {
  stopPlayback();
  state.frameRequestId += 1;
  if (state.frameAbortController) state.frameAbortController.abort();
  state.frameAbortController = null;
  state.playbackLoadPending = false;
  state.run = name;
  state.frameIndex = 0;
  state.conservation = null;
  state.metricIndex = 0;
  ui.diagnosticsState.textContent = "COMPUTING SPHERICAL INTEGRALS…";
  ui.tabs.replaceChildren();
  const metadata = await fetchJSON(`/api/runs/${encodeURIComponent(name)}/metadata`);
  if (state.run !== name) return;
  state.metadata = metadata;
  state.scaleMax = runMaximumSpeed(metadata) ?? 0;
  ui.points.textContent = Number(metadata.grid.point_count).toLocaleString("ja-JP");
  ui.timeline.max = String(metadata.available_steps.length - 1);
  ui.datasetState.textContent = `${metadata.available_frame_count} FRAMES READY`;
  setupViews(metadata);
  await loadFrame(0);
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
    view.renderer.render(view.scene, view.camera);
    views.map.renderer.render(views.map.scene, views.map.camera);
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
  window.addEventListener("resize", drawChart);
  await selectRun(state.runs[0].name);
}

requestAnimationFrame(animate);
init().catch(showFatal);
