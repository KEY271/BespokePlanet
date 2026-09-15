const TWO_PI = Math.PI * 2;

export function wrapLongitude(longitude) {
  return ((longitude + Math.PI) % TWO_PI + TWO_PI) % TWO_PI - Math.PI;
}

function lowerBound(values, target) {
  let low = 0;
  let high = values.length;
  while (low < high) {
    const middle = (low + high) >> 1;
    if (values[middle] < target) low = middle + 1;
    else high = middle;
  }
  return low;
}

/** Interpolate eastward/northward wind on the reduced Gaussian grid. */
export function createWindSampler(grid, eastward, northward) {
  const mu = Float64Array.from(grid.mu, Number);
  const nlon = Uint32Array.from(grid.nlon, Number);
  const offsets = Uint32Array.from(grid.ring_offsets, Number);
  if (eastward.length !== grid.point_count || northward.length !== grid.point_count) {
    throw new Error("Wind component point count does not match metadata");
  }

  const sampleRing = (values, ring, longitude) => {
    const count = nlon[ring];
    const gridCoordinate = ((wrapLongitude(longitude) + Math.PI) / TWO_PI) * count;
    const west = Math.floor(gridCoordinate) % count;
    const fraction = gridCoordinate - Math.floor(gridCoordinate);
    const east = (west + 1) % count;
    const start = offsets[ring];
    return values[start + west] * (1 - fraction) + values[start + east] * fraction;
  };

  return (longitude, latitude) => {
    const targetMu = Math.sin(latitude);
    if (targetMu < mu[0] || targetMu > mu[mu.length - 1]) return null;
    const upper = Math.min(mu.length - 1, Math.max(1, lowerBound(mu, targetMu)));
    const lower = upper - 1;
    const span = mu[upper] - mu[lower];
    const amount = span > 0 ? (targetMu - mu[lower]) / span : 0;
    const u = sampleRing(eastward, lower, longitude) * (1 - amount)
      + sampleRing(eastward, upper, longitude) * amount;
    const v = sampleRing(northward, lower, longitude) * (1 - amount)
      + sampleRing(northward, upper, longitude) * amount;
    const speed = Math.hypot(u, v);
    return Number.isFinite(speed) ? { u, v, speed } : null;
  };
}

function directionAt(sampleWind, longitude, latitude, minimumSpeed) {
  const wind = sampleWind(longitude, latitude);
  const cosLatitude = Math.cos(latitude);
  if (!wind || wind.speed < minimumSpeed || Math.abs(cosLatitude) < 0.04) return null;
  return {
    longitude: wind.u / (wind.speed * cosLatitude),
    latitude: wind.v / wind.speed,
    speed: wind.speed,
  };
}

function integrateHalf(sampleWind, seed, sign, options) {
  const points = [];
  let longitude = seed.longitude;
  let latitude = seed.latitude;
  for (let stepIndex = 0; stepIndex < options.maxSteps; stepIndex += 1) {
    const first = directionAt(sampleWind, longitude, latitude, options.minimumSpeed);
    if (!first) break;
    const signedStep = options.stepRadians * sign;
    const middleLongitude = longitude + first.longitude * signedStep * 0.5;
    const middleLatitude = latitude + first.latitude * signedStep * 0.5;
    if (middleLatitude <= options.minimumLatitude || middleLatitude >= options.maximumLatitude) break;
    const middle = directionAt(sampleWind, middleLongitude, middleLatitude, options.minimumSpeed);
    if (!middle) break;
    const nextLongitude = longitude + middle.longitude * signedStep;
    const nextLatitude = latitude + middle.latitude * signedStep;
    if (nextLatitude <= options.minimumLatitude || nextLatitude >= options.maximumLatitude) break;
    const nextWind = sampleWind(nextLongitude, nextLatitude);
    if (!nextWind || nextWind.speed < options.minimumSpeed) break;
    longitude = nextLongitude;
    latitude = nextLatitude;
    points.push({ longitude, latitude, speed: nextWind.speed });

    if (stepIndex > 24) {
      const latitudeDistance = latitude - seed.latitude;
      const longitudeDistance = wrapLongitude(longitude - seed.longitude) * Math.cos(latitude);
      if (Math.hypot(longitudeDistance, latitudeDistance) < options.stepRadians * 1.4) break;
    }
  }
  return points;
}

export function traceStreamline(sampleWind, seed, options) {
  const seedWind = sampleWind(seed.longitude, seed.latitude);
  if (!seedWind || seedWind.speed < options.minimumSpeed) return [];
  const backward = integrateHalf(sampleWind, seed, -1, options).reverse();
  const forward = integrateHalf(sampleWind, seed, 1, options);
  return [...backward, { ...seed, speed: seedWind.speed }, ...forward];
}

const DENSITY_SETTINGS = {
  low: { latitudeCount: 8, longitudeCount: 12, maxSteps: 64 },
  medium: { latitudeCount: 12, longitudeCount: 18, maxSteps: 78 },
  high: { latitudeCount: 16, longitudeCount: 24, maxSteps: 92 },
};

/** Build instantaneous streamlines. Longitudes stay unwrapped to retain direction at the map seam. */
export function buildStreamlines(grid, eastward, northward, density = "medium") {
  const settings = DENSITY_SETTINGS[density] ?? DENSITY_SETTINGS.medium;
  const sampleWind = createWindSampler(grid, eastward, northward);
  let maximumSpeed = 0;
  for (let index = 0; index < eastward.length; index += 1) {
    maximumSpeed = Math.max(maximumSpeed, Math.hypot(eastward[index], northward[index]));
  }
  if (!(maximumSpeed > 0)) return { lines: [], maximumSpeed: 0 };

  const minimumLatitude = Math.asin(Number(grid.mu[0])) + 1e-5;
  const maximumLatitude = Math.asin(Number(grid.mu.at(-1))) - 1e-5;
  const latitudeMargin = Math.min(0.06, (maximumLatitude - minimumLatitude) * 0.08);
  const options = {
    minimumLatitude,
    maximumLatitude,
    minimumSpeed: Math.max(maximumSpeed * 0.005, 1e-8),
    stepRadians: 0.026,
    maxSteps: settings.maxSteps,
  };
  const lines = [];
  for (let latitudeIndex = 0; latitudeIndex < settings.latitudeCount; latitudeIndex += 1) {
    const fraction = settings.latitudeCount === 1 ? 0.5 : latitudeIndex / (settings.latitudeCount - 1);
    const latitude = minimumLatitude + latitudeMargin
      + fraction * (maximumLatitude - minimumLatitude - 2 * latitudeMargin);
    const offset = latitudeIndex % 2 === 0 ? 0 : 0.5;
    for (let longitudeIndex = 0; longitudeIndex < settings.longitudeCount; longitudeIndex += 1) {
      const longitude = -Math.PI + TWO_PI * (longitudeIndex + offset) / settings.longitudeCount;
      const line = traceStreamline(sampleWind, { longitude, latitude }, options);
      if (line.length >= 8) lines.push(line);
    }
  }
  return { lines, maximumSpeed };
}
