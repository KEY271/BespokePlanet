// Field catalog of the land-sea visualizer: raw output fields with continuous
// colour scales, and the binned/categorical climatology maps of the former R
// analysis (same breaks; palettes are R's hcl.colors values).
import {
  cutIndex, KOPPEN_GROUPS, KOPPEN_TYPES, SNOW_BREAKS, ICE_BREAKS, snowIceClass, snowWaterAtFraction,
  OPEN_OCEAN_LAND_FRACTION,
} from "./climate.js";

// Continuous colour maps, low to high. The sequential ramp runs dark to light so
// that small values recede into the dark surface.
export const COLORMAPS = {
  thermal: ["#002F70", "#2056A2", "#567DCA", "#8DA4DD", "#BBC8ED", "#E1E7F7", "#F9F1F1",
    "#F6CFCF", "#E7A6A6", "#D07979", "#B14949", "#792727", "#5F1415"],
  sequential: ["#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb"],
  diverging: ["#1c5cab", "#3987e5", "#86b6ef", "#383835", "#f0a0a0", "#e66767", "#a83232"],
  terrain: ["#027C1E", "#649334", "#97A753", "#C0B878", "#DEC79D", "#F1D5BF", "#E2E2E2"],
};
export const MISSING_COLOR = "#252c34";

const K = (value) => value - 273.15;

// kind: thermal | sequential | diverging | terrain. mask: land (f_L > 0),
// ocean (f_L < 1) or ice (f_L < 1 and A > 0); masked cells are drawn grey.
export const RAW_FIELDS = {
  surface_temperature: { label: "地表温度", unit: "°C", kind: "thermal", convert: K },
  land_temperature: { label: "陸温度", unit: "°C", kind: "thermal", convert: K, mask: "land" },
  surface_air_temperature: { label: "地上気温", unit: "°C", kind: "thermal", convert: K },
  ocean_temperature: { label: "海洋混合層の水温", unit: "°C", kind: "thermal", convert: K, mask: "ocean" },
  deep_temperature: { label: "陸の深層温度", unit: "°C", kind: "thermal", convert: K, mask: "land" },
  surface_pressure: { label: "地表気圧", unit: "hPa", kind: "sequential", convert: (v) => v / 100 },
  log_surface_pressure: { label: "地表気圧", unit: "hPa", kind: "sequential", convert: (v) => Math.exp(v) / 100 },
  precipitation: { label: "降水量", unit: "mm/日", kind: "sequential", zeroBased: true },
  evaporation: { label: "蒸発量", unit: "mm/日", kind: "sequential", zeroBased: true },
  precipitable_water: { label: "可降水量", unit: "kg/m²", kind: "sequential", zeroBased: true },
  cloud_cover: { label: "雲量", unit: "1", kind: "sequential", zeroBased: true },
  surface_water: { label: "陸のバケツの水", unit: "kg/m²", kind: "sequential", zeroBased: true, mask: "land" },
  surface_wetness: { label: "陸面の湿潤度", unit: "1", kind: "sequential", zeroBased: true, mask: "land" },
  runoff: { label: "陸の流出", unit: "mm/日", kind: "sequential", zeroBased: true, mask: "land" },
  sea_ice_fraction: { label: "海氷面積率", unit: "1", kind: "sequential", zeroBased: true, mask: "ocean" },
  sea_ice_volume: { label: "海氷体積（海面積当たり）", unit: "m", kind: "sequential", zeroBased: true, mask: "ocean" },
  sea_ice_thickness: { label: "海氷の厚さ", unit: "m", kind: "sequential", zeroBased: true, mask: "ice" },
  sea_ice_temperature: { label: "海氷の表面温度", unit: "°C", kind: "thermal", convert: K, mask: "ice" },
  snow_water: { label: "陸の積雪（水当量）", unit: "kg/m²", kind: "sequential", zeroBased: true, mask: "land" },
  snow_fraction: { label: "陸の積雪率", unit: "1", kind: "sequential", zeroBased: true, mask: "land" },
  snowfall: { label: "降雪量", unit: "mm/日", kind: "sequential", zeroBased: true },
  snow_melt: { label: "陸の融雪量", unit: "mm/日", kind: "sequential", zeroBased: true, mask: "land" },
  temperature: { label: "気温", unit: "°C", kind: "thermal", convert: K },
  specific_humidity: { label: "比湿", unit: "g/kg", kind: "sequential", convert: (v) => v * 1000 },
  u: { label: "東西風 u", unit: "m/s", kind: "diverging" },
  v: { label: "南北風 v", unit: "m/s", kind: "diverging" },
  speed: { label: "風速", unit: "m/s", kind: "sequential", zeroBased: true },
  zeta: { label: "相対渦度 ζ", unit: "s⁻¹", kind: "diverging" },
  delta: { label: "発散 δ", unit: "s⁻¹", kind: "diverging" },
  land_fraction: { label: "陸面率 f_L", unit: "1", kind: "sequential", zeroBased: true },
  surface_height: { label: "地表高度 z_s", unit: "m", kind: "terrain" },
  ocean_q_flux: { label: "海洋の Q フラックス", unit: "W/m²", kind: "diverging" },
};

export const MONTHLY_ORDER = [
  "surface_temperature", "surface_air_temperature", "land_temperature", "ocean_temperature", "deep_temperature", "surface_pressure",
  "precipitation", "evaporation", "precipitable_water", "cloud_cover", "surface_water", "surface_wetness", "runoff",
  "snowfall", "snow_water", "snow_fraction", "snow_melt",
  "sea_ice_fraction", "sea_ice_volume", "sea_ice_thickness", "sea_ice_temperature",
];
export const YEARLY_SURFACE_ORDER = [
  "surface_temperature", "land_temperature", "ocean_temperature", "deep_temperature", "log_surface_pressure",
  "cloud_cover", "surface_water", "snow_water", "snow_fraction",
  "sea_ice_fraction", "sea_ice_volume", "sea_ice_thickness", "sea_ice_temperature",
];
export const YEARLY_LEVEL_ORDER = ["temperature", "u", "v", "speed", "specific_humidity", "zeta", "delta"];
export const STATIC_ORDER = ["land_fraction", "surface_height", "ocean_q_flux"];

const OCEAN = "#b9d9eb";
const LAND = "#e6e2d8";
const entries = (colors, labels) => colors.map((color, i) => ({ color, label: labels[i] }));
const binned = (value, breaks, right, clamp = null) => {
  const x = clamp ? Math.min(clamp[1], Math.max(clamp[0], value)) : value;
  return cutIndex(x, breaks, right);
};

export const KOPPEN_GROUP_COLORS = { A: "#2f9e44", B: "#d8b365", C: "#ffd43b", D: "#4dabf7", E: "#f1f3f5" };
export const KOPPEN_GROUP_LABELS = { A: "A 熱帯", B: "B 乾燥帯", C: "C 温帯", D: "D 冷帯", E: "E 寒帯" };
export const KOPPEN_TYPE_COLORS = {
  Af: "#0000fe", Am: "#0078ff", Aw: "#46aafa",
  BWh: "#ff0000", BWk: "#ff9695", BSh: "#f5a301", BSk: "#ffdb63",
  Csa: "#ffff00", Csb: "#c6c700", Csc: "#969600",
  Cwa: "#96ff96", Cwb: "#63c764", Cwc: "#329633",
  Cfa: "#c6ff4e", Cfb: "#66ff33", Cfc: "#33c701",
  Dsa: "#ff00fe", Dsb: "#c600c7", Dsc: "#963295", Dsd: "#966495",
  Dwa: "#abb1ff", Dwb: "#5a77db", Dwc: "#4c51b5", Dwd: "#320087",
  Dfa: "#00ffff", Dfb: "#37c8ff", Dfc: "#007e7d", Dfd: "#00456e",
  ET: "#b2b2b2", EF: "#686868",
};

const SNOW_COLORS = ["#9B92DA", "#8175CB", "#6759B2", "#4C3F8F", "#312271"];
const ICE_COLORS = ["#659FDA", "#2C86CA", "#006BAC", "#005089", "#00366C"];
const percentLabels = (breaks) => breaks.slice(0, -1).map((b, i) => `${Math.round(100 * b)}–${Math.min(100, Math.round(100 * breaks[i + 1]))}%`);

function snowIceLegend(maskingWaterEquivalent) {
  const format = (x) => (!Number.isFinite(x) ? "" : x < 10 ? x.toFixed(1) : x.toFixed(0));
  const snow = percentLabels(SNOW_BREAKS).map((label, i) => {
    const lower = format(snowWaterAtFraction(SNOW_BREAKS[i], maskingWaterEquivalent));
    const upper = format(snowWaterAtFraction(SNOW_BREAKS[i + 1], maskingWaterEquivalent));
    return `積雪率 ${label}（S ${upper ? `${lower}–${upper}` : `≥ ${lower}`} kg/m²）`;
  });
  return [
    { color: LAND, label: "陸・積雪率 1% 未満" },
    ...entries(SNOW_COLORS, snow),
    { color: "#dde6ec", label: "水面・海氷 1% 未満" },
    ...entries(ICE_COLORS, percentLabels(ICE_BREAKS).map((label) => `海氷面積率 ${label}`)),
  ];
}

const monthColumn = (climate, month) => climate.calendarMonths.indexOf(month);

/**
 * Climatology maps. legend: entries {color,label}; classify(climate, i) returns
 * the legend entry of cell i. value(climate, i) is shown in the hover readout.
 * The land threshold enters through climate.land.
 */
export function climatologyFields(context) {
  const iceFraction = [-0.001, 0.001, 0.15, 0.3, 0.5, 0.7, 0.85, 1.001];
  const iceFractionColors = [OCEAN, "#F9F9F9", "#D0E4FF", "#99BFEF", "#5295D4", "#0066A5", "#00366C"];
  const iceFractionMap = (month) => ({
    label: `海氷面積率（${month}月）`,
    unit: "1",
    legend: [...entries(iceFractionColors, ["0", "0–15%", "15–30%", "30–50%", "50–70%", "70–85%", "85–100%"]), { color: LAND, label: "陸" }],
    classify: (c, i) => {
      if (c.land[i]) return 7;
      const column = monthColumn(c, month);
      return binned(c.seaIceFractionMonthly[column][i], iceFraction, true, [0, 1]);
    },
    value: (c, i) => c.seaIceFractionMonthly[monthColumn(c, month)][i],
  });
  const snowIceMap = (month) => ({
    label: `積雪と海氷（${month}月）`,
    unit: "1",
    legend: snowIceLegend(context.maskingWaterEquivalent),
    classify: (c, i) => {
      const column = monthColumn(c, month);
      return snowIceClass(c.land[i], c.snowFraction[column][i], c.seaIceFractionMonthly[column][i]);
    },
    value: (c, i) => {
      const column = monthColumn(c, month);
      return c.land[i] ? c.snowFraction[column][i] : c.seaIceFractionMonthly[column][i];
    },
    valueLabel: (c, i) => (c.land[i] ? "積雪率 f" : "海氷面積率 A"),
  });
  const landOnly = (spec) => ({
    ...spec,
    legend: [...spec.legend, { color: OCEAN, label: "海" }],
    classify: (c, i) => (c.land[i] ? spec.classify(c, i) : spec.legend.length),
    value: (c, i) => (c.land[i] ? spec.value?.(c, i) : undefined),
  });
  const presentTypes = (c) => KOPPEN_TYPES.filter((_, t) => {
    for (let i = 0; i < c.koppenType.length; i += 1) if (c.land[i] && c.koppenType[i] === t) return true;
    return false;
  });

  return {
    koppen_type: landOnly({
      label: "ケッペン–ガイガーの気候型（Peel et al. 2007、C/D 境界 −3 °C、地上気温）",
      legend: KOPPEN_TYPES.map((code) => ({ color: KOPPEN_TYPE_COLORS[code], label: code })),
      legendFilter: (c) => new Set([...presentTypes(c), "海"]),
      classify: (c, i) => c.koppenType[i],
      value: (c, i) => KOPPEN_TYPES[c.koppenType[i]],
    }),
    koppen_group: landOnly({
      label: "ケッペンの気候群（C/D 境界 −3 °C、地上気温）",
      legend: KOPPEN_GROUPS.map((g) => ({ color: KOPPEN_GROUP_COLORS[g], label: KOPPEN_GROUP_LABELS[g] })),
      classify: (c, i) => c.koppenGroup[i],
      value: (c, i) => KOPPEN_GROUPS[c.koppenGroup[i]],
    }),
    terrain: landOnly({
      label: "陸の標高",
      unit: "m",
      legend: entries(["#027C1E", "#649334", "#97A753", "#C0B878", "#DEC79D", "#F1D5BF", "#E2E2E2"],
        ["250 m 未満", "250–500 m", "500–1000 m", "1000–1500 m", "1500–2000 m", "2000–2500 m", "2500 m 以上"]),
      classify: (c, i) => binned(c.surfaceHeight[i], [-Infinity, 250, 500, 1000, 1500, 2000, 2500, Infinity], false),
      value: (c, i) => c.surfaceHeight[i],
    }),
    land_temperature_annual: landOnly({
      label: "陸の年平均地表温度",
      unit: "°C",
      legend: entries(["#002F70", "#3C6FC4", "#9BAEE2", "#DFE5F7", "#F9DFDF", "#E39D9D", "#B74F4F", "#5F1415"],
        ["−30 °C 未満", "−30〜−20 °C", "−20〜−10 °C", "−10〜0 °C", "0〜10 °C", "10〜20 °C", "20〜30 °C", "30 °C 以上"]),
      classify: (c, i) => binned(c.landTemperatureC[i], [-Infinity, -30, -20, -10, 0, 10, 20, 30, Infinity], false),
      value: (c, i) => c.landTemperatureC[i],
    }),
    air_temperature_annual: landOnly({
      label: "陸の年平均地上気温",
      unit: "°C",
      legend: entries(["#002F70", "#3C6FC4", "#9BAEE2", "#DFE5F7", "#F9DFDF", "#E39D9D", "#B74F4F", "#5F1415"],
        ["−30 °C 未満", "−30〜−20 °C", "−20〜−10 °C", "−10〜0 °C", "0〜10 °C", "10〜20 °C", "20〜30 °C", "30 °C 以上"]),
      classify: (c, i) => binned(c.airTemperatureC[i], [-Infinity, -30, -20, -10, 0, 10, 20, 30, Infinity], false),
      value: (c, i) => c.airTemperatureC[i],
    }),
    precipitation_annual: landOnly({
      label: "陸の年降水量",
      unit: "mm/年",
      legend: entries(["#FCFFDD", "#DEF6D2", "#AEE4C1", "#65CCB3", "#00ADAE", "#0088B1", "#005396", "#26185F"],
        ["250 mm 未満", "250–500 mm", "500–1000 mm", "1000–1500 mm", "1500–2000 mm", "2000–3000 mm", "3000–4000 mm", "4000 mm 以上"]),
      classify: (c, i) => binned(c.precipitation[i], [-Infinity, 250, 500, 1000, 1500, 2000, 3000, 4000, Infinity], false),
      value: (c, i) => c.precipitation[i],
    }),
    surface_water_annual: landOnly({
      label: "陸のバケツの水の年平均",
      unit: "kg/m²",
      legend: entries(["#FCFFDD", "#E6F9D5", "#C6EDC8", "#9ADCBB", "#59C8B2", "#00B1AE", "#0095AF", "#0073AD", "#00468B", "#26185F"],
        Array.from({ length: 10 }, (_, i) => `${15 * i}–${15 * (i + 1)} kg/m²`)),
      classify: (c, i) => binned(c.surfaceWater[i], Array.from({ length: 11 }, (_, i) => 15 * i), true, [0, 150]),
      value: (c, i) => c.surfaceWater[i],
    }),
    p_minus_e_annual: landOnly({
      label: "陸の年間 P − E（陸面が捨てる流出）",
      unit: "mm/年",
      legend: entries(["#533600", "#986D00", "#D1A358", "#F2D4B0", "#F6F6F6", "#B8DFD9", "#67B2A9", "#007B71", "#003D34"],
        ["−500 未満", "−500〜−250", "−250〜−100", "−100〜0", "0〜100", "100〜250", "250〜500", "500〜1000", "1000 以上"]),
      classify: (c, i) => binned(c.pMinusE[i], [-Infinity, -500, -250, -100, 0, 100, 250, 500, 1000, Infinity], false),
      value: (c, i) => c.pMinusE[i],
    }),
    cloud_cover_annual: {
      label: "年平均の実効雲量",
      unit: "1",
      legend: entries(["#F9F9F9", "#E4F0FF", "#CBE0FF", "#ADCCF6", "#8CB6E9", "#659FDA", "#2C86CA", "#006BAC", "#005089", "#00366C"],
        Array.from({ length: 10 }, (_, i) => `${10 * i}–${10 * (i + 1)}%`)),
      classify: (c, i) => binned(c.cloudCover[i], Array.from({ length: 11 }, (_, i) => i / 10), true, [0, 1]),
      value: (c, i) => c.cloudCover[i],
    },
    precipitation_season: {
      label: "降水量の JJA − DJF",
      unit: "mm/日",
      legend: entries(["#5F1415", "#9D3D3D", "#CB6F70", "#E6A4A4", "#F7D3D3", "#F6F6F6", "#D3DBF4", "#A3B4E5", "#6889D0", "#265BAB", "#002F70"],
        ["−8 未満", "−8〜−4", "−4〜−2", "−2〜−1", "−1〜−0.25", "−0.25〜0.25", "0.25〜1", "1〜2", "2〜4", "4〜8", "8 以上"].map((s) => `${s} mm/日`)),
      classify: (c, i) => binned(c.precipitationSeason[i], [-Infinity, -8, -4, -2, -1, -0.25, 0.25, 1, 2, 4, 8, Infinity], false),
      value: (c, i) => c.precipitationSeason[i],
    },
    pressure_season: {
      label: "地表気圧の JJA − DJF",
      unit: "hPa",
      legend: entries(["#492050", "#82498C", "#B574C2", "#D2A9DB", "#E8D4ED", "#F1F1F1", "#C8E1C9", "#91C392", "#4E9D4F", "#256C26", "#023903"],
        ["−12 未満", "−12〜−8", "−8〜−4", "−4〜−2", "−2〜−0.5", "−0.5〜0.5", "0.5〜2", "2〜4", "4〜8", "8〜12", "12 以上"].map((s) => `${s} hPa`)),
      classify: (c, i) => binned(c.pressureSeason[i], [-Infinity, -12, -8, -4, -2, -0.5, 0.5, 2, 4, 8, 12, Infinity], false),
      value: (c, i) => c.pressureSeason[i],
    },
    ocean_pressure_anomaly: {
      label: `外洋の地表気圧偏差（f_L ≤ ${OPEN_OCEAN_LAND_FRACTION}）`,
      unit: "hPa",
      legend: [...entries(["#002F70", "#265BAB", "#6889D0", "#A3B4E5", "#D3DBF4", "#F6F6F6", "#F7D3D3", "#E6A4A4", "#CB6F70", "#9D3D3D", "#5F1415"],
        ["−12 未満", "−12〜−8", "−8〜−4", "−4〜−2", "−2〜−1", "−1〜1", "1〜2", "2〜4", "4〜8", "8〜12", "12 以上"].map((s) => `${s} hPa`)),
      { color: "#f2f2f2", label: "外洋以外" }],
      classify: (c, i) => (c.openOcean[i]
        ? binned(c.oceanPressureAnomaly[i], [-Infinity, -12, -8, -4, -2, -1, 1, 2, 4, 8, 12, Infinity], false) : 11),
      value: (c, i) => c.oceanPressureAnomaly[i],
    },
    sea_ice_march: iceFractionMap(3),
    sea_ice_september: iceFractionMap(9),
    sea_ice_thickness_annual: {
      label: "海氷の厚さ（氷面積で重み付けした年平均）",
      unit: "m",
      legend: [...entries(["#FFFFC8", "#FCE7A0", "#F7C252", "#F39300", "#EB5500", "#BE1D00", "#7D0025"],
        ["0.25 m 未満", "0.25–0.5 m", "0.5–1 m", "1–2 m", "2–4 m", "4–8 m", "8 m 以上"]),
      { color: OCEAN, label: "海氷なし" }, { color: LAND, label: "陸" }],
      classify: (c, i) => {
        if (c.land[i]) return 8;
        if (!(c.iceFraction[i] > 0)) return 7;
        return binned(c.iceThickness[i], [-0.001, 0.25, 0.5, 1, 2, 4, 8, Infinity], true);
      },
      value: (c, i) => c.iceThickness[i],
    },
    sea_ice_temperature_annual: {
      label: "海氷の表面温度（氷面積で重み付けした年平均）",
      unit: "°C",
      legend: [...entries(["#002F70", "#517AC9", "#B4C2EB", "#F6F6F6", "#EDB4B5", "#C05D5D", "#5F1415"],
        ["−40 °C 未満", "−40〜−30 °C", "−30〜−20 °C", "−20〜−10 °C", "−10〜−5 °C", "−5〜0 °C", "0 °C 以上"]),
      { color: OCEAN, label: "海氷なし" }, { color: LAND, label: "陸" }],
      classify: (c, i) => {
        if (c.land[i]) return 8;
        if (!(c.iceFraction[i] > 0) || !Number.isFinite(c.iceTemperatureC[i])) return 7;
        return binned(c.iceTemperatureC[i], [-Infinity, -40, -30, -20, -10, -5, 0, Infinity], true);
      },
      value: (c, i) => c.iceTemperatureC[i],
    },
    snow_ice_march: snowIceMap(3),
    snow_ice_september: snowIceMap(9),
  };
}

export const CLIMATOLOGY_ORDER = [
  "koppen_type", "koppen_group", "terrain", "air_temperature_annual", "land_temperature_annual", "precipitation_annual",
  "surface_water_annual", "p_minus_e_annual", "cloud_cover_annual", "precipitation_season", "pressure_season",
  "ocean_pressure_anomaly", "snow_ice_march", "snow_ice_september", "sea_ice_march", "sea_ice_september",
  "sea_ice_thickness_annual", "sea_ice_temperature_annual",
];

export function hexToRgb(hex) {
  const value = Number.parseInt(hex.slice(1), 16);
  return [((value >> 16) & 255) / 255, ((value >> 8) & 255) / 255, (value & 255) / 255];
}

/** Colour of t in [0, 1] on a continuous colour map, as linear-free sRGB triplet. */
export function sampleColormap(stops, t) {
  const x = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const lower = Math.min(stops.length - 2, Math.floor(x));
  const amount = x - lower;
  const a = hexToRgb(stops[lower]);
  const b = hexToRgb(stops[lower + 1]);
  return [a[0] + (b[0] - a[0]) * amount, a[1] + (b[1] - a[1]) * amount, a[2] + (b[2] - a[2]) * amount];
}

export function cssGradient(stops) {
  return `linear-gradient(90deg, ${stops.join(", ")})`;
}
