import assert from "node:assert/strict";
import test from "node:test";

import * as C from "./static/climate.js";

const metadata = { grid: { mu: [-1 / Math.sqrt(3), 1 / Math.sqrt(3)], nlon: [4, 4] } };
const grid = C.createGrid(metadata);
const close = (actual, expected, tolerance = 1e-9) => assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}`);

test("grid follows the solver convention: k = 0 at longitude 0, south to north", () => {
  assert.equal(grid.pointCount, 8);
  assert.deepEqual(Array.from(grid.lon.slice(0, 4)), [0, 90, 180, 270]);
  assert.deepEqual(Array.from(grid.mapLon.slice(0, 4)), [0, 90, -180, -90]);
  assert.ok(grid.lat[0] < 0 && grid.lat[4] > 0);
  close(grid.weight.reduce((a, b) => a + b, 0), 1);
  assert.equal(C.cellAt(grid, 91, 20), 5);
  assert.equal(C.cellAt(grid, -89, -20), 3);
});

test("Gauss-Legendre weights of the nodes integrate polynomials", () => {
  const mu = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
  const w = C.gaussLegendreWeights(mu);
  close(w.reduce((a, b) => a + b, 0), 2, 1e-12);
  close(w.reduce((sum, weight, i) => sum + weight * mu[i] ** 2, 0), 2 / 3, 1e-12);
});

test("calendar months and complete years", () => {
  assert.equal(C.calendarMonthOf(1, 4), 4);
  assert.equal(C.calendarMonthOf(10, 4), 1);
  assert.equal(C.startMonthOf({ simulation: { start_calendar_time: "0001-04-01 00:00:00" } }), 4);
  const years = C.completeYears([...Array.from({ length: 25 }, (_, i) => i + 1)], 12);
  assert.deepEqual(years.map((year) => year.year), [1, 2]);
  assert.deepEqual(years[1].months.slice(0, 2), [13, 14]);
});

test("cutIndex matches R's cut with include.lowest", () => {
  const breaks = [-Infinity, 0, 10, Infinity];
  assert.equal(C.cutIndex(-5, breaks), 0);
  assert.equal(C.cutIndex(0, breaks), 1);
  assert.equal(C.cutIndex(10, breaks), 2);
  const closed = [0, 15, 30];
  assert.equal(C.cutIndex(0, closed, true), 0);
  assert.equal(C.cutIndex(15, closed, true), 0);
  assert.equal(C.cutIndex(15.1, closed, true), 1);
  assert.equal(C.cutIndex(30, closed, true), 1);
  assert.equal(C.cutIndex(Number.NaN, closed, true), -1);
  assert.equal(C.findInterval(0.5, [0.01, 0.1, 0.3, 0.5]), 4);
});

test("coastline of a single land cell is closed and inside the map", () => {
  const land = new Uint8Array(8);
  land[5] = 1;
  const segments = C.coastlineSegments(grid, land);
  const length = segments.reduce((sum, [x0, y0, x1, y1]) => sum + Math.hypot(x1 - x0, y1 - y0), 0);
  // Two meridional edges (the ring spans the equator to 90N) and the southern
  // edge; the pole is not an edge.
  close(length, 2 * 90 + 90);
  for (const [x0, , x1] of segments) assert.ok(Math.abs(x0) <= 180 && Math.abs(x1) <= 180);
  land.fill(0);
  land[6] = 1; // centred on the seam: split into both map edges
  const seam = C.coastlineSegments(grid, land);
  assert.ok(seam.every(([x0, , x1]) => Math.abs(x0) <= 180 && Math.abs(x1) <= 180));
  close(seam.reduce((sum, [x0, y0, x1, y1]) => sum + Math.hypot(x1 - x0, y1 - y0), 0), 270);
});

test("site snapping uses great-circle distance to land cells", () => {
  const land = Uint8Array.from([0, 1, 0, 0, 0, 0, 0, 1]);
  const [site] = C.snapSites(grid, [{ site: "x", lon: 100, lat: -10 }], land);
  assert.equal(site.index, 1);
  assert.ok(site.snapDistance > 0);
});

function uniformInput({ temperatureC, precipitation }) {
  const n = grid.pointCount;
  const months = 12;
  const constant = (value) => Array.from({ length: months }, (_, m) => new Float64Array(n).fill(typeof value === "function" ? value(m) : value));
  return {
    landFraction: Float64Array.from([1, 1, 1, 1, 0, 0, 0, 0.6]),
    surfaceHeight: new Float64Array(n),
    monthly: {
      surface_temperature: constant((m) => temperatureC(m) + 273.15),
      land_temperature: constant((m) => temperatureC(m) + 273.15),
      ocean_temperature: constant(290),
      precipitation: constant(precipitation),
      evaporation: constant(1),
      cloud_cover: constant(0.5),
      surface_pressure: constant(100000),
      surface_water: constant(75),
      sea_ice_fraction: constant(0.5),
      sea_ice_volume: constant(1),
      sea_ice_temperature: constant(260),
      snow_fraction: constant(0),
    },
    zonalV: Array.from({ length: months }, () => Float64Array.from([1, 1, 1, 1])),
  };
}

const options = (threshold = 0.5) => ({
  threshold,
  calendarMonths: Array.from({ length: 12 }, (_, m) => C.calendarMonthOf(m + 1, 4)),
  outputMonths: Array.from({ length: 12 }, (_, m) => m + 1),
  daysPerMonth: 30,
  aHalf: [0, 100, 0],
  bHalf: [0, 0.5, 1],
  referenceHalfPressure: [0, 50000, 100000],
  gravity: 9.80616,
  terrain: "analytic",
});

test("Köppen groups: wet tropics, dry desert, polar", () => {
  const tropical = C.computeClimate(grid, uniformInput({ temperatureC: () => 26, precipitation: 10 }), options());
  assert.equal(C.KOPPEN_GROUPS[tropical.koppenGroup[0]], "A");
  assert.equal(C.KOPPEN_TYPES[tropical.koppenType[0]], "Af");
  const desert = C.computeClimate(grid, uniformInput({ temperatureC: () => 26, precipitation: 0.1 }), options());
  assert.equal(C.KOPPEN_TYPES[desert.koppenType[0]], "BWh");
  const polar = C.computeClimate(grid, uniformInput({ temperatureC: () => -20, precipitation: 1 }), options());
  assert.equal(C.KOPPEN_TYPES[polar.koppenType[0]], "EF");
  // A cold month of -2 degC is C with the -3 degC boundary and D with the 0 degC one.
  const boundary = C.computeClimate(grid, uniformInput({ temperatureC: (m) => (m === 9 ? -2 : 15), precipitation: 3 }), options());
  assert.equal(C.KOPPEN_GROUPS[boundary.koppenGroup[0]], "C");
  assert.equal(C.KOPPEN_GROUPS[boundary.koppenGroupZero[0]], "D");
});

test("land threshold decides which cells count as land", () => {
  const input = uniformInput({ temperatureC: () => 26, precipitation: 10 });
  assert.equal(C.computeClimate(grid, input, options(0.5)).land[7], 1);
  const strict = C.computeClimate(grid, input, options(0.7));
  assert.equal(strict.land[7], 0);
  const threshold = strict.summary.find((row) => row.metric === "land_threshold");
  assert.equal(threshold.value, 0.7);
});

test("mass streamfunction integrates v dp from the top", () => {
  const input = uniformInput({ temperatureC: () => 10, precipitation: 1 });
  const climate = C.computeClimate(grid, input, options());
  const factor = 2 * Math.PI * C.EARTH_RADIUS_M / 9.80616 * grid.cosLatitude[0];
  const dp = [100 + 0.5 * 100000, -100 + 0.5 * 100000];
  close(climate.streamfunction.psi[0][0], factor * 0.5 * dp[0], 1e-3);
  close(climate.streamfunction.psi[0][1], factor * (dp[0] + 0.5 * dp[1]), 1e-3);
});

test("sea-ice integrals use ocean-area weights", () => {
  const climate = C.computeClimate(grid, uniformInput({ temperatureC: () => 10, precipitation: 1 }), options());
  const oceanFraction = grid.weight.reduce((sum, w, i) => sum + w * (1 - [1, 1, 1, 1, 0, 0, 0, 0.6][i]), 0);
  const global = climate.seaIceMonthly.find((row) => row.region === "global");
  close(global.areaM2, 0.5 * oceanFraction * 4 * Math.PI * C.EARTH_RADIUS_M ** 2, 1);
  close(global.thicknessM, 2);
  assert.ok(Number.isNaN(climate.iceThickness[0]));
  close(climate.iceTemperatureC[4], 260 - 273.15);
});

test("snow and sea-ice classes", () => {
  assert.equal(C.snowIceClass(true, 0.005, 0), 0);
  assert.equal(C.snowIceClass(true, 1, 0), 5);
  assert.equal(C.snowIceClass(false, 0, 0.005), 6);
  assert.equal(C.snowIceClass(false, 0, 0.9), 11);
  close(C.snowWaterAtFraction(0.5, 50), 50);
});
