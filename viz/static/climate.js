// Climate diagnostics of the land-sea output, ported from the former R analysis
// scripts (analyze_moist_land_sea*.R, snow_sea_ice_map.R). Pure functions without
// DOM or WebGL access so that they can be tested with `node --test`.
//
// Grid convention (docs/dynamics/octahedral-gaussian-grid.md): rings are ordered
// south to north, point k of ring j sits at longitude 360 k / nlon_j degrees, and
// flat fields are ring-major with longitude varying fastest.

export const EARTH_RADIUS_M = 6.371e6; // core/src/config/planet_parameters.f90
export const OPEN_OCEAN_LAND_FRACTION = 0.01;
export const MONTH_ABBREVIATIONS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
export const MONTH_NAMES = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
export const KOPPEN_GROUP_NAMES = { A: "熱帯", B: "乾燥帯", C: "温帯", D: "冷帯", E: "寒帯" };

export const KOPPEN_GROUPS = ["A", "B", "C", "D", "E"];
// C/D boundary on the coldest-month temperature, used by both the groups and the types.
export const KOPPEN_COLD_BOUNDARY_C = -3;
export const KOPPEN_TYPES = [
  "Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk",
  "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc",
  "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd",
  "ET", "EF",
];

// Reference sites. Each is snapped to the nearest cell whose land fraction is
// at least the land threshold.
export const ANALYTIC_SITES = [
  { site: "A north interior", name: "大陸 A 北部内陸", lon: 45, lat: 50 },
  { site: "A south interior", name: "大陸 A 南部内陸", lon: 45, lat: 30 },
  { site: "B equatorial interior", name: "大陸 B 赤道内陸", lon: 280, lat: -5 },
  { site: "B south interior", name: "大陸 B 南部内陸", lon: 280, lat: -30 },
  { site: "C interior", name: "大陸 C 内陸", lon: 150, lat: -30 },
  { site: "South polar cap", name: "南極冠", lon: 180, lat: -80 },
];
// The Earth Köppen type is the observed type at the real site, for reference only.
export const EARTH_SITES = [
  { site: "Manaus (Amazon)", name: "マナウス（アマゾン）", lon: 300.0, lat: -3.1, earth: "Af" },
  { site: "Kinshasa (Congo)", name: "キンシャサ（コンゴ）", lon: 15.3, lat: -4.3, earth: "Aw" },
  { site: "Kolkata (India)", name: "コルカタ（インド）", lon: 88.4, lat: 22.6, earth: "Aw" },
  { site: "Cairo (Sahara)", name: "カイロ（サハラ）", lon: 31.2, lat: 30.0, earth: "BWh" },
  { site: "Alice Springs (Australia)", name: "アリススプリングス（オーストラリア）", lon: 133.9, lat: -23.7, earth: "BWh" },
  { site: "Buenos Aires (Pampas)", name: "ブエノスアイレス（パンパ）", lon: 301.6, lat: -34.6, earth: "Cfa" },
  { site: "Paris (W Europe)", name: "パリ（西ヨーロッパ）", lon: 2.3, lat: 48.9, earth: "Cfb" },
  { site: "Chicago (N America)", name: "シカゴ（北アメリカ）", lon: 272.4, lat: 41.9, earth: "Dfa" },
  { site: "Yakutsk (Siberia)", name: "ヤクーツク（シベリア）", lon: 129.7, lat: 62.0, earth: "Dfd" },
  { site: "Lhasa (Tibet)", name: "ラサ（チベット）", lon: 91.1, lat: 29.7, earth: "Dwb" },
  { site: "Greenland interior", name: "グリーンランド内陸", lon: 318.0, lat: 72.0, earth: "EF" },
  { site: "Vostok (Antarctica)", name: "ボストーク（南極）", lon: 106.8, lat: -78.5, earth: "EF" },
];
// Land boxes of the Earth-terrain checks (longitude 0..360 east).
export const EARTH_BOXES = {
  south_asia: [70, 100, 5, 30],
  east_asia: [100, 125, 20, 45],
  amazon: [290, 310, -15, 5],
  tibet_upwind: [80, 95, 20, 28],
  tibet_lee: [80, 95, 36, 45],
  andes_west: [283, 289, -35, -20],
  andes_east: [293, 302, -35, -20],
};

const DEG = Math.PI / 180;

/** Gauss-Legendre weights recovered from the nodes, as in the R scripts. */
export function gaussLegendreWeights(mu) {
  const n = mu.length;
  return mu.map((x) => {
    let previous = 1;
    let current = x;
    for (let k = 2; k <= n; k += 1) {
      const next = ((2 * k - 1) * x * current - (k - 1) * previous) / k;
      previous = current;
      current = next;
    }
    if (n === 1) previous = 1;
    const derivative = n * (x * current - previous) / (x * x - 1);
    return 2 / ((1 - x * x) * derivative * derivative);
  });
}

export function createGrid(metadata) {
  const mu = metadata.grid.mu.map(Number);
  const nlon = metadata.grid.nlon.map(Number);
  const nlat = mu.length;
  const offsets = [0];
  for (const count of nlon) offsets.push(offsets.at(-1) + count);
  const pointCount = offsets.at(-1);
  const latitude = mu.map((value) => Math.asin(value) / DEG);
  const latEdges = [-90];
  for (let j = 1; j < nlat; j += 1) latEdges.push(0.5 * (latitude[j - 1] + latitude[j]));
  latEdges.push(90);
  const gaussian = gaussLegendreWeights(mu);
  const ring = new Int32Array(pointCount);
  const lon = new Float64Array(pointCount);
  const mapLon = new Float64Array(pointCount);
  const lat = new Float64Array(pointCount);
  const weight = new Float64Array(pointCount);
  for (let j = 0; j < nlat; j += 1) {
    for (let k = 0; k < nlon[j]; k += 1) {
      const index = offsets[j] + k;
      ring[index] = j;
      lon[index] = (360 * k) / nlon[j];
      mapLon[index] = lon[index] >= 180 ? lon[index] - 360 : lon[index];
      lat[index] = latitude[j];
      weight[index] = gaussian[j] / 2 / nlon[j];
    }
  }
  return {
    nlat, nlon, offsets, pointCount, mu, latitude, latEdges,
    cosLatitude: mu.map((value) => Math.sqrt(Math.max(0, 1 - value * value))),
    ring, lon, mapLon, lat, weight,
  };
}

/** Index of the cell containing (longitude, latitude) in degrees, or -1. */
export function cellAt(grid, longitude, latitude) {
  if (!(latitude >= -90 && latitude <= 90) || !Number.isFinite(longitude)) return -1;
  let j = 0;
  while (j < grid.nlat - 1 && latitude > grid.latEdges[j + 1]) j += 1;
  const count = grid.nlon[j];
  const east = ((longitude % 360) + 360) % 360;
  const k = Math.round((east / 360) * count) % count;
  return grid.offsets[j] + k;
}

/** Calendar month (1-12) of output month m (m0001 starts at the start date). */
export function calendarMonthOf(month, startMonth = 4) {
  return ((startMonth - 1 + month - 1) % 12 + 12) % 12 + 1;
}

export function startMonthOf(metadata) {
  const match = /^\d+-(\d{2})-/.exec(String(metadata.simulation?.start_calendar_time ?? ""));
  return match ? Number(match[1]) : 4;
}

/** Complete simulation years: year y holds output months (y-1)P+1 .. yP. */
export function completeYears(months, monthsPerYear) {
  const available = new Set(months);
  const years = [];
  const last = Math.max(0, ...months);
  for (let year = 1; year * monthsPerYear <= last; year += 1) {
    const indices = Array.from({ length: monthsPerYear }, (_, i) => (year - 1) * monthsPerYear + i + 1);
    if (indices.every((month) => available.has(month))) years.push({ year, months: indices });
  }
  return years;
}

/**
 * Class of `value` for the breaks of R's cut(..., include.lowest = TRUE).
 * right = false: [b0,b1), ..., [b(n-1),bn]; right = true: [b0,b1], ..., (b(n-1),bn].
 * Returns -1 for NaN or values outside the breaks.
 */
export function cutIndex(value, breaks, right = false) {
  if (!Number.isFinite(value) && !(value === Infinity || value === -Infinity)) return -1;
  const last = breaks.length - 1;
  if (value < breaks[0] || value > breaks[last]) return -1;
  if (right) {
    if (value <= breaks[1]) return 0;
    for (let i = 1; i < last; i += 1) if (value <= breaks[i + 1]) return i;
    return last - 1;
  }
  if (value >= breaks[last - 1]) return last - 1;
  for (let i = 0; i < last; i += 1) if (value < breaks[i + 1]) return i;
  return last - 1;
}

/** R's findInterval: the number of breaks that are <= value. */
export function findInterval(value, breaks) {
  let count = 0;
  while (count < breaks.length && value >= breaks[count]) count += 1;
  return count;
}

export function landMask(landFraction, threshold) {
  const land = new Uint8Array(landFraction.length);
  for (let i = 0; i < land.length; i += 1) land[i] = landFraction[i] >= threshold ? 1 : 0;
  return land;
}

export function ringMean(grid, values) {
  const result = new Float64Array(grid.nlat);
  for (let j = 0; j < grid.nlat; j += 1) {
    let sum = 0;
    for (let index = grid.offsets[j]; index < grid.offsets[j + 1]; index += 1) sum += values[index];
    result[j] = sum / grid.nlon[j];
  }
  return result;
}

function meanOf(arrays) {
  const result = new Float64Array(arrays[0].length);
  for (const values of arrays) for (let i = 0; i < result.length; i += 1) result[i] += values[i];
  for (let i = 0; i < result.length; i += 1) result[i] /= arrays.length;
  return result;
}

/** Annual mean of ring means, i.e. rowMeans of the monthly zonal means. */
export function zonalAnnualMean(grid, monthly) {
  return meanOf(monthly.map((values) => ringMean(grid, values)));
}

export function weightedMean(values, weights, include = null) {
  let sum = 0;
  let total = 0;
  for (let i = 0; i < values.length; i += 1) {
    if (include && !include[i]) continue;
    sum += weights[i] * values[i];
    total += weights[i];
  }
  return total > 0 ? sum / total : Number.NaN;
}

/** Great-circle snap of target sites to the nearest candidate cell. */
export function snapSites(grid, targets, candidate) {
  return targets.map((target) => {
    const targetLon = target.lon * DEG;
    const targetLat = target.lat * DEG;
    let best = -1;
    let bestCosine = -Infinity;
    for (let i = 0; i < grid.pointCount; i += 1) {
      if (!candidate[i]) continue;
      const cosine = Math.sin(targetLat) * Math.sin(grid.lat[i] * DEG)
        + Math.cos(targetLat) * Math.cos(grid.lat[i] * DEG) * Math.cos(grid.lon[i] * DEG - targetLon);
      if (cosine > bestCosine) { bestCosine = cosine; best = i; }
    }
    return {
      ...target,
      index: best,
      snapDistance: best < 0 ? Number.NaN : Math.acos(Math.min(1, bestCosine)) / DEG,
    };
  });
}

/**
 * Coastline of a land mask: the edges between land and non-land cells, as
 * segments [lon0, lat0, lon1, lat1] in map longitude (-180..180) degrees.
 */
export function coastlineSegments(grid, land) {
  const segments = [];
  const cellIndex = (j, i) => grid.offsets[j] + ((i % grid.nlon[j]) + grid.nlon[j]) % grid.nlon[j];
  for (let j = 0; j < grid.nlat; j += 1) {
    const width = 360 / grid.nlon[j];
    const south = grid.latEdges[j];
    const north = grid.latEdges[j + 1];
    for (let i = 0; i < grid.nlon[j]; i += 1) {
      const index = cellIndex(j, i);
      if (!land[index]) continue;
      const left = grid.mapLon[index] - width / 2;
      const right = grid.mapLon[index] + width / 2;
      if (!land[cellIndex(j, i + 1)]) segments.push([right, south, right, north]);
      if (!land[cellIndex(j, i - 1)]) segments.push([left, south, left, north]);
      for (const neighbour of [j - 1, j + 1]) {
        if (neighbour < 0 || neighbour >= grid.nlat) continue;
        const neighbourWidth = 360 / grid.nlon[neighbour];
        // Every neighbouring-ring cell overlapping this cell's longitude span.
        const first = Math.floor((grid.lon[index] - width / 2) / neighbourWidth + 0.5);
        const last = Math.ceil((grid.lon[index] + width / 2) / neighbourWidth - 0.5);
        const edge = neighbour < j ? south : north;
        for (let n = first; n <= last; n += 1) {
          const other = cellIndex(neighbour, n);
          if (land[other]) continue;
          const otherCentre = n * neighbourWidth;
          const from = Math.max(grid.lon[index] - width / 2, otherCentre - neighbourWidth / 2);
          const to = Math.min(grid.lon[index] + width / 2, otherCentre + neighbourWidth / 2);
          if (to - from <= 1e-9) continue;
          const offset = grid.mapLon[index] - grid.lon[index];
          segments.push([from + offset, edge, to + offset, edge]);
        }
      }
    }
  }
  return segments.flatMap(wrapSegment);
}

/** Keep segments inside -180..180, splitting those that cross the seam. */
function wrapSegment([x0, y0, x1, y1]) {
  const shift = (x) => (x < -180 ? x + 360 : x > 180 ? x - 360 : x);
  if (x0 >= -180 && x0 <= 180 && x1 >= -180 && x1 <= 180) return [[x0, y0, x1, y1]];
  if (x0 === x1) return [[shift(x0), y0, shift(x1), y1]];
  const [left, right] = x0 < x1 ? [x0, x1] : [x1, x0];
  if (left < -180) return [[left + 360, y0, 180, y1], [-180, y0, right, y1]].filter(([a, , b]) => b - a > 1e-9);
  return [[left, y0, 180, y1], [-180, y0, right - 360, y1]].filter(([a, , b]) => b - a > 1e-9);
}

/**
 * Köppen inputs shared by the group and type classifications. The temperature is
 * the monthly near-surface air temperature (K), not the land skin temperature.
 */
function koppenInputs(grid, airTemperatureK, precipitationMmDay, calendarMonths, daysPerMonth) {
  const n = grid.pointCount;
  const months = calendarMonths.length;
  const out = {
    annualTemperature: new Float64Array(n), warmest: new Float64Array(n), coldest: new Float64Array(n),
    monthsAbove10: new Int32Array(n), annualPrecipitation: new Float64Array(n), summerFraction: new Float64Array(n),
    driestMonth: new Float64Array(n), driestSummer: new Float64Array(n), wettestSummer: new Float64Array(n),
    driestWinter: new Float64Array(n), wettestWinter: new Float64Array(n), dryness: new Float64Array(n),
  };
  const northSummer = calendarMonths.map((month) => month >= 4 && month <= 9);
  for (let i = 0; i < n; i += 1) {
    let sumT = 0; let warmest = -Infinity; let coldest = Infinity; let above10 = 0;
    let annualP = 0; let summerP = 0; let driest = Infinity;
    let driestSummer = Infinity; let wettestSummer = -Infinity; let driestWinter = Infinity; let wettestWinter = -Infinity;
    const north = grid.lat[i] >= 0;
    for (let m = 0; m < months; m += 1) {
      const t = airTemperatureK[m][i] - 273.15;
      const p = precipitationMmDay[m][i] * daysPerMonth;
      sumT += t;
      warmest = Math.max(warmest, t);
      coldest = Math.min(coldest, t);
      if (t > 10) above10 += 1;
      annualP += p;
      driest = Math.min(driest, p);
      if (north === northSummer[m]) {
        summerP += p;
        driestSummer = Math.min(driestSummer, p);
        wettestSummer = Math.max(wettestSummer, p);
      } else {
        driestWinter = Math.min(driestWinter, p);
        wettestWinter = Math.max(wettestWinter, p);
      }
    }
    const annualT = sumT / months;
    const fraction = annualP > 0 ? summerP / annualP : 0.5;
    const adjustment = fraction >= 0.7 ? 280 : fraction <= 0.3 ? 0 : 140;
    out.annualTemperature[i] = annualT;
    out.warmest[i] = warmest;
    out.coldest[i] = coldest;
    out.monthsAbove10[i] = above10;
    out.annualPrecipitation[i] = annualP;
    out.summerFraction[i] = fraction;
    out.driestMonth[i] = driest;
    out.driestSummer[i] = driestSummer;
    out.wettestSummer[i] = wettestSummer;
    out.driestWinter[i] = driestWinter;
    out.wettestWinter[i] = wettestWinter;
    out.dryness[i] = 20 * annualT + adjustment;
  }
  return out;
}

function koppenGroupAt(inputs, i) {
  if (inputs.annualPrecipitation[i] < inputs.dryness[i]) return 1; // B
  if (inputs.coldest[i] >= 18) return 0; // A
  if (inputs.warmest[i] < 10) return 4; // E
  if (inputs.coldest[i] > KOPPEN_COLD_BOUNDARY_C) return 2; // C
  return 3; // D
}

/** Full Köppen-Geiger type of Peel et al. (2007), with the -3 degC C/D boundary of the groups. */
function koppenTypeAt(inputs, i, groupIndex) {
  const group = KOPPEN_GROUPS[groupIndex];
  const annualP = inputs.annualPrecipitation[i];
  let code = group;
  if (group === "B") {
    code = `B${annualP < 0.5 * inputs.dryness[i] ? "W" : "S"}${inputs.annualTemperature[i] >= 18 ? "h" : "k"}`;
  } else if (group === "A") {
    if (inputs.driestMonth[i] >= 60) code = "Af";
    else if (inputs.driestMonth[i] >= 100 - annualP / 25) code = "Am";
    else code = "Aw";
  } else if (group === "C" || group === "D") {
    let season = "f";
    if (inputs.driestSummer[i] < 40 && inputs.driestSummer[i] < inputs.wettestWinter[i] / 3) season = "s";
    else if (inputs.driestWinter[i] < inputs.wettestSummer[i] / 10) season = "w";
    let heat = "c";
    if (inputs.warmest[i] >= 22) heat = "a";
    else if (inputs.monthsAbove10[i] >= 4) heat = "b";
    else if (group === "D" && inputs.coldest[i] < -38) heat = "d";
    code = `${group}${season}${heat}`;
  } else if (group === "E") {
    code = inputs.warmest[i] > 0 ? "ET" : "EF";
  }
  return KOPPEN_TYPES.indexOf(code);
}

/** Annual-mean meridional mass streamfunction from monthly zonal-mean v and p_s. */
export function massStreamfunction(grid, zonalV, surfacePressure, aHalf, bHalf, gravity, radius = EARTH_RADIUS_M) {
  const nlev = aHalf.length - 1;
  const psi = Array.from({ length: grid.nlat }, () => new Float64Array(nlev));
  const residual = Array.from({ length: grid.nlat }, () => new Float64Array(zonalV.length));
  zonalV.forEach((v, m) => {
    const zonalPressure = ringMean(grid, surfacePressure[m]);
    for (let j = 0; j < grid.nlat; j += 1) {
      let cumulative = 0;
      let column = 0;
      for (let k = 0; k < nlev; k += 1) {
        const dp = (aHalf[k + 1] - aHalf[k]) + zonalPressure[j] * (bHalf[k + 1] - bHalf[k]);
        const transport = (2 * Math.PI * radius / gravity) * grid.cosLatitude[j] * v[k * grid.nlat + j] * dp;
        cumulative += transport;
        column += transport;
        psi[j][k] += (cumulative - 0.5 * transport) / zonalV.length;
      }
      residual[j][m] = column;
    }
  });
  return { psi, residual };
}

/**
 * Final-year climate of one analysis year.
 *
 * input: { landFraction, surfaceHeight, monthly: {field: Float64Array[months]},
 *          zonalV: Float64Array[months] } with months in output order.
 * options: { threshold, calendarMonths, daysPerMonth, aHalf, bHalf,
 *            referenceHalfPressure, gravity, maskingWaterEquivalent, terrain }
 */
export function computeClimate(grid, input, options) {
  const n = grid.pointCount;
  const { monthly, landFraction } = input;
  const months = options.calendarMonths.length;
  const dpm = options.daysPerMonth;
  const land = landMask(landFraction, options.threshold);
  const mean = (field) => meanOf(monthly[field]);
  const sumTimes = (field, factor) => {
    const out = new Float64Array(n);
    for (const values of monthly[field]) for (let i = 0; i < n; i += 1) out[i] += values[i] * factor;
    return out;
  };
  const toCelsius = (values) => values.map((value) => value - 273.15);

  const surfaceTemperatureC = toCelsius(mean("surface_temperature"));
  const landTemperatureC = toCelsius(mean("land_temperature"));
  const airTemperatureC = toCelsius(mean("surface_air_temperature"));
  const oceanTemperatureC = toCelsius(mean("ocean_temperature"));
  const precipitation = sumTimes("precipitation", dpm);
  const evaporation = sumTimes("evaporation", dpm);
  const cloudCover = mean("cloud_cover");
  const surfacePressureHpa = mean("surface_pressure").map((value) => value / 100);
  const surfaceWater = mean("surface_water");
  const pMinusE = precipitation.map((value, i) => value - evaporation[i]);

  const inputs = koppenInputs(grid, monthly.surface_air_temperature, monthly.precipitation, options.calendarMonths, dpm);
  const koppenGroup = new Int8Array(n);
  const koppenType = new Int8Array(n);
  for (let i = 0; i < n; i += 1) {
    koppenGroup[i] = koppenGroupAt(inputs, i);
    koppenType[i] = koppenTypeAt(inputs, i, koppenGroup[i]);
  }

  // Seasonal contrast (JJA minus DJF).
  const seasonMean = (field, calendar) => {
    const columns = options.calendarMonths.map((month, m) => (calendar.includes(month) ? m : -1)).filter((m) => m >= 0);
    const out = new Float64Array(n);
    for (const m of columns) for (let i = 0; i < n; i += 1) out[i] += monthly[field][m][i] / columns.length;
    return out;
  };
  const jja = [6, 7, 8];
  const djf = [12, 1, 2];
  const precipitationJja = seasonMean("precipitation", jja);
  const precipitationDjf = seasonMean("precipitation", djf);
  const pressureJja = seasonMean("surface_pressure", jja);
  const pressureDjf = seasonMean("surface_pressure", djf);
  const precipitationSeason = precipitationJja.map((value, i) => value - precipitationDjf[i]);
  const pressureSeason = pressureJja.map((value, i) => (value - pressureDjf[i]) / 100);

  // Ocean surface pressure: departure from the open-ocean mean and grid-scale ripple.
  const openOcean = landFraction.map((value) => (value <= OPEN_OCEAN_LAND_FRACTION ? 1 : 0));
  const oceanMeanPressure = weightedMean(surfacePressureHpa, grid.weight, openOcean);
  const oceanPressureAnomaly = surfacePressureHpa.map((value, i) => (openOcean[i] ? value - oceanMeanPressure : Number.NaN));
  const gridScalePressure = new Float64Array(n);
  for (let j = 0; j < grid.nlat; j += 1) {
    const start = grid.offsets[j];
    const count = grid.nlon[j];
    for (let k = 0; k < count; k += 1) {
      const left = surfacePressureHpa[start + (k - 1 + count) % count];
      const right = surfacePressureHpa[start + (k + 1) % count];
      gridScalePressure[start + k] = surfacePressureHpa[start + k] - 0.5 * (left + right);
    }
  }

  // Sea ice: concentration and volume per ocean area.
  const oceanWeight = grid.weight.map((weight, i) => weight * (1 - landFraction[i]));
  const earthArea = 4 * Math.PI * EARTH_RADIUS_M ** 2;
  const iceFraction = mean("sea_ice_fraction");
  const iceVolume = mean("sea_ice_volume");
  const iceThickness = new Float64Array(n);
  const iceTemperatureC = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    let weighted = 0;
    let area = 0;
    for (let m = 0; m < months; m += 1) {
      weighted += monthly.sea_ice_fraction[m][i] * monthly.sea_ice_temperature[m][i];
      area += monthly.sea_ice_fraction[m][i];
    }
    const ocean = landFraction[i] < 1;
    iceThickness[i] = ocean && iceFraction[i] > 0 ? iceVolume[i] / iceFraction[i] : Number.NaN;
    iceTemperatureC[i] = ocean && iceFraction[i] > 0 && area > 0 ? weighted / area - 273.15 : Number.NaN;
  }
  const regions = {
    global: () => true,
    north: (i) => grid.lat[i] >= 0,
    south: (i) => grid.lat[i] < 0,
  };
  const seaIceMonthly = [];
  for (const [region, inRegion] of Object.entries(regions)) {
    for (let m = 0; m < months; m += 1) {
      let area = 0;
      let volume = 0;
      for (let i = 0; i < n; i += 1) {
        if (!inRegion(i)) continue;
        area += oceanWeight[i] * monthly.sea_ice_fraction[m][i];
        volume += oceanWeight[i] * monthly.sea_ice_volume[m][i];
      }
      seaIceMonthly.push({
        region,
        outputMonth: options.outputMonths?.[m] ?? m + 1,
        calendarMonth: options.calendarMonths[m],
        areaM2: area * earthArea,
        volumeM3: volume * earthArea,
        thicknessM: area > 0 ? volume / area : Number.NaN,
      });
    }
  }

  // Zonal means over land and ocean.
  const zonal = {
    latitude: grid.latitude,
    surfaceTemperatureK: zonalAnnualMean(grid, monthly.surface_temperature),
    precipitation: zonalAnnualMean(grid, monthly.precipitation),
    evaporation: zonalAnnualMean(grid, monthly.evaporation),
    cloudCover: zonalAnnualMean(grid, monthly.cloud_cover),
  };

  const streamfunction = massStreamfunction(
    grid, input.zonalV, monthly.surface_pressure, options.aHalf, options.bHalf, options.gravity,
  );
  const half = options.referenceHalfPressure;
  const referencePressureHpa = half.slice(0, -1).map((value, k) => 0.5 * (value + half[k + 1]) / 100);

  const climate = {
    land, threshold: options.threshold, calendarMonths: options.calendarMonths,
    outputMonths: options.outputMonths ?? [],
    surfaceTemperatureC, landTemperatureC, airTemperatureC, oceanTemperatureC, precipitation, evaporation, pMinusE,
    cloudCover, surfacePressureHpa, surfaceWater, koppenGroup, koppenType,
    koppenInputs: inputs, precipitationSeason, pressureSeason, openOcean, oceanPressureAnomaly,
    gridScalePressure, iceFraction, iceVolume, iceThickness, iceTemperatureC, oceanWeight,
    seaIceMonthly, zonal, streamfunction, referencePressureHpa,
    surfaceHeight: input.surfaceHeight, landFraction,
    monthlyAirTemperature: monthly.surface_air_temperature, monthlyPrecipitation: monthly.precipitation,
    snowFraction: monthly.snow_fraction, seaIceFractionMonthly: monthly.sea_ice_fraction,
  };
  climate.sites = snapSites(grid, options.terrain === "earth" ? EARTH_SITES : ANALYTIC_SITES, land);
  climate.summary = summarize(grid, climate, options);
  return climate;
}

const maxOf = (values) => values.reduce((a, b) => (b > a ? b : a), -Infinity);
const minOf = (values) => values.reduce((a, b) => (b < a ? b : a), Infinity);

function summarize(grid, climate, options) {
  const n = grid.pointCount;
  const { land, landFraction } = climate;
  const w = grid.weight;
  const landWeight = w.map((weight, i) => weight * landFraction[i]);
  const oceanWeight = w.map((weight, i) => weight * (1 - landFraction[i]));
  const psiValues = climate.streamfunction.psi.flatMap((row) => Array.from(row, (value) => value / 1e9));
  const residual = climate.streamfunction.residual.map((row) => Math.abs(row.reduce((a, b) => a + b, 0) / row.length));
  const landIndices = [];
  for (let i = 0; i < n; i += 1) if (land[i]) landIndices.push(i);
  const landArea = landIndices.reduce((sum, i) => sum + landWeight[i], 0);
  const rows = [];
  const add = (metric, value, label) => rows.push({ metric, value, label });
  const months = climate.outputMonths;
  if (months.length) {
    add("final_month_start", months[0], "最初の出力月");
    add("final_month_end", months.at(-1), "最後の出力月");
  }
  add("land_threshold", climate.threshold, "陸の閾値 f_L");
  add("global_land_area_fraction", landWeight.reduce((a, b) => a + b, 0), "全球平均の陸面率");
  add("land_cell_area_fraction", landIndices.reduce((sum, i) => sum + w[i], 0), "陸とみなしたセルの面積割合");
  add("maximum_streamfunction_1e9_kg_s", maxOf(psiValues), "質量流線関数の最大（10⁹ kg/s）");
  add("minimum_streamfunction_1e9_kg_s", minOf(psiValues), "質量流線関数の最小（10⁹ kg/s）");
  add("maximum_surface_mass_residual_1e9_kg_s", maxOf(residual) / 1e9, "鉛直積分した質量輸送の残差の最大（10⁹ kg/s）");
  add("global_area_mean_cloud_cover", weightedMean(climate.cloudCover, w), "雲量（全球）");
  add("land_area_mean_cloud_cover", weightedMean(climate.cloudCover, landWeight), "雲量（陸面積平均）");
  add("ocean_area_mean_cloud_cover", weightedMean(climate.cloudCover, oceanWeight), "雲量（海面積平均）");
  add("global_area_mean_surface_temperature_c", weightedMean(climate.surfaceTemperatureC, w), "地表温度（全球、°C）");
  add("land_area_mean_surface_temperature_c", weightedMean(climate.landTemperatureC, landWeight), "陸温度（陸面積平均、°C）");
  add("land_area_mean_surface_air_temperature_c", weightedMean(climate.airTemperatureC, landWeight), "地上気温（陸面積平均、°C）");
  add("ocean_area_mean_surface_temperature_c", weightedMean(climate.oceanTemperatureC, oceanWeight), "海水温（海面積平均、°C）");
  add("land_area_mean_surface_water_kg_m-2", weightedMean(climate.surfaceWater, landWeight), "陸のバケツの水（kg/m²）");
  add("land_precipitation_minus_evaporation_mm_yr", weightedMean(climate.pMinusE, landWeight), "P − E（陸面積平均、mm/年）");
  add("ocean_precipitation_minus_evaporation_mm_yr", weightedMean(climate.pMinusE, oceanWeight), "P − E（海面積平均、mm/年）");
  const landT = landIndices.map((i) => climate.landTemperatureC[i]);
  add("minimum_land_annual_temperature_c", landT.length ? minOf(landT) : Number.NaN, "最も寒い陸セル（°C）");
  add("maximum_land_annual_temperature_c", landT.length ? maxOf(landT) : Number.NaN, "最も暑い陸セル（°C）");
  const polar = (test) => weightedMean(climate.landTemperatureC, w, land.map((value, i) => value && test(grid.lat[i])));
  add("antarctic_land_annual_temperature_c", polar((lat) => lat < -65), "南緯65°以南の陸（°C）");
  add("arctic_land_annual_temperature_c", polar((lat) => lat > 65), "北緯65°以北の陸（°C）");
  const polarAir = (test) => weightedMean(climate.airTemperatureC, w, land.map((value, i) => value && test(grid.lat[i])));
  add("antarctic_land_annual_surface_air_temperature_c", polarAir((lat) => lat < -65), "南緯65°以南の陸の地上気温（°C）");
  add("arctic_land_annual_surface_air_temperature_c", polarAir((lat) => lat > 65), "北緯65°以北の陸の地上気温（°C）");
  const anomaly = [];
  const anomalyWeight = [];
  const ripple = [];
  for (let i = 0; i < n; i += 1) {
    if (!climate.openOcean[i]) continue;
    anomaly.push(climate.oceanPressureAnomaly[i]);
    anomalyWeight.push(w[i]);
    ripple.push(Math.abs(climate.gridScalePressure[i]));
  }
  const anomalyMean = weightedMean(anomaly, anomalyWeight);
  add("ocean_surface_pressure_anomaly_max_abs_hpa", maxOf(anomaly.map(Math.abs)), "外洋の地表気圧偏差の最大絶対値（hPa）");
  add("ocean_surface_pressure_anomaly_sd_hpa",
    Math.sqrt(weightedMean(anomaly.map((value) => (value - anomalyMean) ** 2), anomalyWeight)), "外洋の地表気圧偏差の標準偏差（hPa）");
  add("ocean_surface_pressure_grid_scale_max_abs_hpa", maxOf(ripple), "外洋の格子スケールの気圧リップルの最大（hPa）");

  if (options.terrain === "earth") {
    const box = ([lonMin, lonMax, latMin, latMax]) => land.map((value, i) => value
      && grid.lat[i] >= latMin && grid.lat[i] <= latMax && grid.lon[i] >= lonMin && grid.lon[i] <= lonMax);
    const boxPrecipitation = (mask, calendar) => {
      const columns = climate.calendarMonths.map((month, m) => (calendar === null || calendar.includes(month) ? m : -1)).filter((m) => m >= 0);
      const values = new Float64Array(n);
      for (const m of columns) for (let i = 0; i < n; i += 1) values[i] += climate.monthlyPrecipitation[m][i] * 360 / columns.length;
      return weightedMean(values, w, mask);
    };
    const boxes = Object.fromEntries(Object.entries(EARTH_BOXES).map(([key, value]) => [key, box(value)]));
    add("south_asia_jja_minus_djf_precipitation_mm_yr",
      boxPrecipitation(boxes.south_asia, [6, 7, 8]) - boxPrecipitation(boxes.south_asia, [12, 1, 2]), "南アジアの降水量 JJA − DJF（mm/年）");
    add("east_asia_jja_minus_djf_precipitation_mm_yr",
      boxPrecipitation(boxes.east_asia, [6, 7, 8]) - boxPrecipitation(boxes.east_asia, [12, 1, 2]), "東アジアの降水量 JJA − DJF（mm/年）");
    add("amazon_djf_minus_jja_precipitation_mm_yr",
      boxPrecipitation(boxes.amazon, [12, 1, 2]) - boxPrecipitation(boxes.amazon, [6, 7, 8]), "アマゾンの降水量 DJF − JJA（mm/年）");
    add("tibet_upwind_annual_precipitation_mm_yr", boxPrecipitation(boxes.tibet_upwind, null), "チベット風上側の降水量（mm/年）");
    add("tibet_lee_annual_precipitation_mm_yr", boxPrecipitation(boxes.tibet_lee, null), "チベット風下側の降水量（mm/年）");
    add("andes_west_annual_precipitation_mm_yr", boxPrecipitation(boxes.andes_west, null), "アンデス西側の降水量（mm/年）");
    add("andes_east_annual_precipitation_mm_yr", boxPrecipitation(boxes.andes_east, null), "アンデス東側の降水量（mm/年）");
  }
  KOPPEN_GROUPS.forEach((group, g) => {
    let area = 0;
    for (const i of landIndices) if (climate.koppenGroup[i] === g) area += landWeight[i];
    add(`land_area_fraction_group_${group}`, landArea > 0 ? area / landArea : Number.NaN, `気候群 ${group} の陸面積割合`);
  });
  const global = climate.seaIceMonthly.filter((row) => row.region === "global");
  let iceVolume = 0;
  let iceArea = 0;
  for (let i = 0; i < n; i += 1) {
    iceVolume += climate.oceanWeight[i] * climate.iceVolume[i];
    iceArea += climate.oceanWeight[i] * climate.iceFraction[i];
  }
  add("final_year_mean_sea_ice_area_m2", global.reduce((sum, row) => sum + row.areaM2, 0) / global.length, "海氷面積の年平均（m²）");
  add("final_year_mean_sea_ice_volume_m3", global.reduce((sum, row) => sum + row.volumeM3, 0) / global.length, "海氷体積の年平均（m³）");
  add("final_year_mean_sea_ice_thickness_m", iceArea > 0 ? iceVolume / iceArea : Number.NaN, "海氷の平均の厚さ（m）");
  add("final_year_max_sea_ice_area_m2", maxOf(global.map((row) => row.areaM2)), "月平均海氷面積の最大（m²）");
  return rows;
}

/**
 * Site record for a climograph: January-December monthly near-surface air
 * temperature and precipitation, the same inputs as the Köppen classification.
 */
export function climographOf(climate, index) {
  const order = [];
  for (let month = 1; month <= 12; month += 1) order.push(climate.calendarMonths.indexOf(month));
  const source = climate.monthlyAirTemperature;
  const temperature = order.map((m) => (m < 0 ? Number.NaN : source[m][index] - 273.15));
  const precipitation = order.map((m) => (m < 0 ? Number.NaN : climate.monthlyPrecipitation[m][index] * 30));
  return {
    index,
    temperature,
    precipitation,
    annualTemperature: climate.airTemperatureC[index],
    annualPrecipitation: climate.precipitation[index],
    koppenGroup: KOPPEN_GROUPS[climate.koppenGroup[index]],
    koppenType: KOPPEN_TYPES[climate.koppenType[index]] ?? "—",
  };
}

// Snow and sea-ice map classes (snow_sea_ice_map.R). Cover below 1% counts as none.
export const SNOW_BREAKS = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0001];
export const ICE_BREAKS = [0.01, 0.15, 0.3, 0.5, 0.85, 1.0001];

/** Class 0 = bare land, 1..5 = snow; 6 = open water, 7..11 = sea ice. */
export function snowIceClass(isLand, snowFraction, iceFraction) {
  if (isLand) return Math.min(5, findInterval(snowFraction, SNOW_BREAKS));
  return 6 + Math.min(5, findInterval(iceFraction, ICE_BREAKS));
}

/** Snow water equivalent S = S_0 f / (1 - f) at a snow-cover fraction f. */
export function snowWaterAtFraction(fraction, maskingWaterEquivalent) {
  return fraction >= 1 ? Infinity : maskingWaterEquivalent * fraction / (1 - fraction);
}
