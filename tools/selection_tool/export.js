function sanitizeBaseName(name) {
  const stripped = name.replace(/\.[^.]+$/, "");
  return stripped || "spectrum";
}

function normalizeRegion(region, index) {
  if (!region || typeof region !== "object") {
    throw new Error(`Region ${index + 1} is not an object.`);
  }

  if (region.type !== "inclusion" && region.type !== "exclusion") {
    throw new Error(`Region ${index + 1} has an invalid type.`);
  }

  const start = Number(region.start);
  const end = Number(region.end);
  if (!Number.isFinite(start) || !Number.isFinite(end)) {
    throw new Error(`Region ${index + 1} has non-numeric bounds.`);
  }

  return {
    type: region.type,
    start: Math.min(start, end),
    end: Math.max(start, end)
  };
}

export function buildRegionsDocument(regions, sourceName = null) {
  return {
    format: "selection_tool_regions",
    version: 1,
    sourceName: sourceName || null,
    regions: regions
      .map((region, index) => normalizeRegion(region, index))
      .sort((a, b) => a.start - b.start)
  };
}

export function buildRegionsJson(regions, sourceName = null) {
  return JSON.stringify(buildRegionsDocument(regions, sourceName), null, 2);
}

export function parseRegionsJson(text) {
  let parsed;
  try {
    parsed = JSON.parse(text);
  } catch (error) {
    throw new Error("The selected file is not valid JSON.");
  }

  const rawRegions = Array.isArray(parsed) ? parsed : parsed && parsed.regions;
  if (!Array.isArray(rawRegions)) {
    throw new Error("The JSON file must contain a top-level regions array.");
  }

  return {
    sourceName:
      !Array.isArray(parsed) && parsed && typeof parsed.sourceName === "string"
        ? parsed.sourceName
        : null,
    regions: rawRegions.map((region, index) => normalizeRegion(region, index))
  };
}

export function exportRegions(regions, sourceName) {
  const text = buildRegionsJson(regions, sourceName);
  const blob = new Blob([text], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  const baseName = sanitizeBaseName(sourceName || "regions");

  link.href = url;
  link.download = `${baseName}_regions.json`;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);

  return `Exported ${regions.length} region${regions.length === 1 ? "" : "s"} to ${link.download}.`;
}
