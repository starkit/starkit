import { parseSpectrumText } from "./parser.js";
import { exportRegions, parseRegionsJson } from "./export.js";
import { ATLAS_CATALOG, DEFAULT_ATLAS_ID, getAtlasById } from "./atlases.js";
import {
  clampRegion,
  createRegion,
  hitTestRegions,
  normalizeBounds,
  regionWidth,
  resizeRegion,
  sortRegions
} from "./regions.js";
import { SpectrumPlot } from "./plot.js";

const MIN_REGION_PIXEL_WIDTH = 4;
const MAX_Y_SAMPLE_POINTS = 2000;
const AUTO_Y_LOWER_QUANTILE = 0.01;
const AUTO_Y_UPPER_QUANTILE = 0.99;
const Y_SCALE_LABELS = {
  auto: "Auto",
  zero: "Zero",
  full: "Full"
};
const NO_ATLAS_ID = "none";
const defaultAtlas = getAtlasById(DEFAULT_ATLAS_ID);

const state = {
  spectrum: null,
  regions: [],
  selectedRegionId: null,
  activeToolType: "inclusion",
  interaction: {
    mode: null,
    regionId: null,
    anchorWavelength: null,
    startWavelength: null,
    originalStart: null,
    originalEnd: null
  },
  draftRegion: null,
  view: null,
  status: "Ready.",
  warnings: [],
  hoverWavelength: null,
  regionCounter: 1,
  zoomMode: false,
  zoomBox: null,
  selectedAtlasId: defaultAtlas.id,
  activeAtlas: defaultAtlas,
  wavelengthUnit: "angstrom",
  regionSourceName: null,
  yScaleMode: "auto"
};

const elements = {
  fileInput: document.querySelector("#spectrum-file"),
  regionsInput: document.querySelector("#regions-file"),
  fileName: document.querySelector("#file-name"),
  dropZone: document.querySelector("#drop-zone"),
  plotCanvas: document.querySelector("#plot-canvas"),
  emptyState: document.querySelector("#empty-state"),
  cursorReadout: document.querySelector("#cursor-readout"),
  summaryPoints: document.querySelector("#summary-points"),
  summaryRange: document.querySelector("#summary-range"),
  summaryFlux: document.querySelector("#summary-flux"),
  summaryRegions: document.querySelector("#summary-regions"),
  yScaleAuto: document.querySelector("#y-scale-auto"),
  yScaleZero: document.querySelector("#y-scale-zero"),
  yScaleFull: document.querySelector("#y-scale-full"),
  unitAngstrom: document.querySelector("#unit-angstrom"),
  unitNm: document.querySelector("#unit-nm"),
  unitRaw: document.querySelector("#unit-raw"),
  statusBox: document.querySelector("#status-box"),
  regionsBody: document.querySelector("#regions-body"),
  selectionLabel: document.querySelector("#selection-label"),
  modeInclusion: document.querySelector("#mode-inclusion"),
  modeExclusion: document.querySelector("#mode-exclusion"),
  deleteSelected: document.querySelector("#delete-selected"),
  zoomBoxToggle: document.querySelector("#zoom-box-toggle"),
  zoomIn: document.querySelector("#zoom-in"),
  zoomOut: document.querySelector("#zoom-out"),
  zoomSelected: document.querySelector("#zoom-selected"),
  atomicAtlasSelect: document.querySelector("#atomic-atlas-select"),
  resetView: document.querySelector("#reset-view"),
  exportRegions: document.querySelector("#export-regions")
};

const plot = new SpectrumPlot(elements.plotCanvas);

function formatNumber(value) {
  if (!Number.isFinite(value)) {
    return "-";
  }
  if (Math.abs(value) >= 1000 || (Math.abs(value) > 0 && Math.abs(value) < 0.01)) {
    return value.toExponential(4);
  }
  return value.toFixed(6).replace(/\.?0+$/, "");
}

function formatPlainDecimal(value) {
  if (!Number.isFinite(value)) {
    return "-";
  }
  return value.toFixed(3).replace(/\.?0+$/, "");
}

function convertWavelength(value) {
  if (!Number.isFinite(value)) {
    return value;
  }
  if (state.wavelengthUnit === "nm") {
    return value / 10;
  }
  return value;
}

function displayToRawWavelength(value) {
  if (!Number.isFinite(value)) {
    return value;
  }
  if (state.wavelengthUnit === "nm") {
    return value * 10;
  }
  return value;
}

function formatWavelength(value) {
  const converted = convertWavelength(value);
  if (state.wavelengthUnit === "angstrom" || state.wavelengthUnit === "nm") {
    return formatPlainDecimal(converted);
  }
  return formatNumber(converted);
}

function getWavelengthUnitLabel() {
  if (state.wavelengthUnit === "angstrom") {
    return "Å";
  }
  if (state.wavelengthUnit === "nm") {
    return "nm";
  }
  return "raw";
}

function describeAtlas(atlas) {
  if (!atlas || !Number.isFinite(atlas.wavelengthMin) || !Number.isFinite(atlas.wavelengthMax)) {
    return "No atomic labels.";
  }

  return `${atlas.label} atlas (${atlas.labels.length} labels, ${formatWavelength(atlas.wavelengthMin)} - ${formatWavelength(atlas.wavelengthMax)} ${getWavelengthUnitLabel()})`;
}

function getAtomicLabelsStatusMessage() {
  if (!state.activeAtlas) {
    return "Atomic labels hidden.";
  }

  return `Atomic labels set to ${state.activeAtlas.label} (${state.activeAtlas.labels.length} labels available).`;
}

function setActiveAtlas(atlasId) {
  if (atlasId === NO_ATLAS_ID) {
    state.selectedAtlasId = NO_ATLAS_ID;
    state.activeAtlas = null;
    return;
  }

  const atlas = getAtlasById(atlasId);
  state.selectedAtlasId = atlas.id;
  state.activeAtlas = atlas;
}

function setStatus(message, extraWarnings = []) {
  state.status = message;
  state.warnings = extraWarnings;
  renderStatus();
}

function renderStatus() {
  const lines = [state.status];
  if (state.warnings.length > 0) {
    lines.push("");
    lines.push(...state.warnings.map((warning) => `Warning: ${warning}`));
  }
  elements.statusBox.textContent = lines.join("\n");
}

function getSpectrumRange() {
  if (!state.spectrum) {
    return null;
  }
  return {
    xMin: state.spectrum.wavelengthMin,
    xMax: state.spectrum.wavelengthMax,
    yMin: state.spectrum.fluxMin,
    yMax: state.spectrum.fluxMax
  };
}

function getVisibleFluxValues(xMin, xMax) {
  if (!state.spectrum) {
    return [];
  }

  const values = [];
  for (const point of state.spectrum.points) {
    if (point.wavelength < xMin) {
      continue;
    }
    if (point.wavelength > xMax) {
      break;
    }
    if (Number.isFinite(point.flux)) {
      values.push(point.flux);
    }
  }

  if (values.length > 0) {
    return values;
  }

  return state.spectrum.points
    .map((point) => point.flux)
    .filter((value) => Number.isFinite(value));
}

function sampleValues(values, maxCount) {
  if (values.length <= maxCount) {
    return values;
  }

  const sampled = [];
  const step = (values.length - 1) / (maxCount - 1);
  for (let index = 0; index < maxCount; index += 1) {
    sampled.push(values[Math.round(index * step)]);
  }
  return sampled;
}

function percentile(sortedValues, fraction) {
  if (sortedValues.length === 0) {
    return 0;
  }
  if (sortedValues.length === 1) {
    return sortedValues[0];
  }

  const position = fraction * (sortedValues.length - 1);
  const lowerIndex = Math.floor(position);
  const upperIndex = Math.ceil(position);
  if (lowerIndex === upperIndex) {
    return sortedValues[lowerIndex];
  }

  const weight = position - lowerIndex;
  return sortedValues[lowerIndex] * (1 - weight) + sortedValues[upperIndex] * weight;
}

function minMax(values) {
  let min = values[0];
  let max = values[0];
  for (let index = 1; index < values.length; index += 1) {
    const value = values[index];
    if (value < min) {
      min = value;
    }
    if (value > max) {
      max = value;
    }
  }
  return { min, max };
}

function padYRange(yMin, yMax) {
  if (!Number.isFinite(yMin) || !Number.isFinite(yMax)) {
    return { yMin: -1, yMax: 1 };
  }

  let lower = yMin;
  let upper = yMax;
  if (lower === upper) {
    const halfSpan = Math.max(Math.abs(lower) * 0.1, 0.1);
    lower -= halfSpan;
    upper += halfSpan;
  }

  const span = upper - lower;
  const padding = Math.max(span * 0.06, Math.abs((upper + lower) / 2) * 0.02, 0.02);
  return {
    yMin: lower - padding,
    yMax: upper + padding
  };
}

function getYRangeForXRange(xMin, xMax, mode) {
  const fluxValues = getVisibleFluxValues(xMin, xMax);
  if (fluxValues.length === 0) {
    return { yMin: -1, yMax: 1 };
  }

  let yMin;
  let yMax;
  if (mode === "full" || fluxValues.length < 8) {
    const range = minMax(fluxValues);
    yMin = range.min;
    yMax = range.max;
  } else {
    const sorted = [...sampleValues(fluxValues, MAX_Y_SAMPLE_POINTS)].sort((a, b) => a - b);
    yMin = percentile(sorted, AUTO_Y_LOWER_QUANTILE);
    yMax = percentile(sorted, AUTO_Y_UPPER_QUANTILE);
  }

  if (mode === "zero") {
    yMin = Math.min(yMin, 0);
    yMax = Math.max(yMax, 0);
  }

  return padYRange(yMin, yMax);
}

function updateViewForXRange(xMin, xMax) {
  if (!state.spectrum) {
    return;
  }

  const bounds = normalizeBounds(
    clampDataX(xMin),
    clampDataX(xMax)
  );
  const yRange = getYRangeForXRange(bounds.start, bounds.end, state.yScaleMode);
  state.view = {
    xMin: bounds.start,
    xMax: bounds.end,
    yMin: yRange.yMin,
    yMax: yRange.yMax
  };
}

function resetView() {
  const range = getSpectrumRange();
  if (!range) {
    return;
  }
  updateViewForXRange(range.xMin, range.xMax);
}

function getMinRegionWidthData() {
  if (!state.view) {
    return 0;
  }
  const bounds = plot.getPlotBounds();
  const dataPerPixel = (state.view.xMax - state.view.xMin) / Math.max(1, bounds.right - bounds.left);
  return MIN_REGION_PIXEL_WIDTH * dataPerPixel;
}

function getSelectedRegion() {
  return state.regions.find((region) => region.id === state.selectedRegionId) || null;
}

function updateModeButtons() {
  elements.modeInclusion.classList.toggle("is-active", state.activeToolType === "inclusion");
  elements.modeExclusion.classList.toggle("is-active", state.activeToolType === "exclusion");
  elements.zoomBoxToggle.classList.toggle("is-active", state.zoomMode);
  elements.yScaleAuto.classList.toggle("is-active", state.yScaleMode === "auto");
  elements.yScaleZero.classList.toggle("is-active", state.yScaleMode === "zero");
  elements.yScaleFull.classList.toggle("is-active", state.yScaleMode === "full");
  elements.unitAngstrom.classList.toggle("is-active", state.wavelengthUnit === "angstrom");
  elements.unitNm.classList.toggle("is-active", state.wavelengthUnit === "nm");
  elements.unitRaw.classList.toggle("is-active", state.wavelengthUnit === "raw");
}

function renderAtlasControls() {
  elements.atomicAtlasSelect.value = state.selectedAtlasId;
  elements.atomicAtlasSelect.title = describeAtlas(state.activeAtlas);
}

function renderSummary() {
  if (!state.spectrum) {
    elements.fileName.textContent = "No spectrum loaded";
    elements.summaryPoints.textContent = "-";
    elements.summaryRange.textContent = "-";
    elements.summaryFlux.textContent = "-";
    elements.summaryRegions.textContent = String(state.regions.length);
    return;
  }

  elements.fileName.textContent = state.spectrum.sourceName;
  elements.summaryPoints.textContent = formatNumber(state.spectrum.pointCount);
  elements.summaryRange.textContent = `${formatWavelength(state.spectrum.wavelengthMin)} - ${formatWavelength(state.spectrum.wavelengthMax)} ${getWavelengthUnitLabel()}`;
  elements.summaryFlux.textContent = `${formatNumber(state.spectrum.fluxMin)} - ${formatNumber(state.spectrum.fluxMax)}`;
  elements.summaryRegions.textContent = String(state.regions.length);
}

function renderSelectionLabel() {
  const region = getSelectedRegion();
  if (!region) {
    elements.selectionLabel.textContent = "No selection";
    return;
  }
  elements.selectionLabel.textContent = `${region.type} ${formatWavelength(region.start)} to ${formatWavelength(region.end)} ${getWavelengthUnitLabel()}`;
}

function renderRegionsTable() {
  const sorted = sortRegions(state.regions);
  if (sorted.length === 0) {
    elements.regionsBody.innerHTML = `
      <tr class="placeholder-row">
        <td colspan="5">No regions defined</td>
      </tr>
    `;
    return;
  }

  elements.regionsBody.innerHTML = sorted
    .map((region) => {
      const isSelected = region.id === state.selectedRegionId;
      const badgeClass = region.type === "inclusion" ? "badge badge-inclusion" : "badge badge-exclusion";
      return `
        <tr data-region-id="${region.id}" class="${isSelected ? "is-selected" : ""}">
          <td><button type="button" class="${badgeClass} row-select">${region.type}</button></td>
          <td><input type="number" step="any" class="boundary-input" data-field="start" value="${convertWavelength(region.start)}" /></td>
          <td><input type="number" step="any" class="boundary-input" data-field="end" value="${convertWavelength(region.end)}" /></td>
          <td>${formatWavelength(regionWidth(region))}</td>
          <td><button type="button" class="row-button row-delete">Delete</button></td>
        </tr>
      `;
    })
    .join("");
}

function renderPlot() {
  plot.draw({
    spectrum: state.spectrum,
    regions: state.regions,
    selectedRegionId: state.selectedRegionId,
    view: state.view,
    draftRegion: state.draftRegion,
    zoomBox: state.zoomBox,
    hoverWavelength: state.hoverWavelength,
    atomicLabels: state.activeAtlas ? state.activeAtlas.labels : [],
    showAtomicLabels: Boolean(state.activeAtlas),
    wavelengthFormatter: formatWavelength,
    wavelengthUnitLabel: getWavelengthUnitLabel()
  });
}

function renderCursorReadout() {
  elements.cursorReadout.textContent =
    state.hoverWavelength === null ? "-" : `${formatWavelength(state.hoverWavelength)} ${getWavelengthUnitLabel()}`;
}

function renderEmptyState() {
  if (!state.spectrum && state.regions.length > 0) {
    elements.emptyState.textContent = "Load a spectrum to visualize imported regions";
  } else {
    elements.emptyState.textContent = "Load a spectrum to begin";
  }
  elements.emptyState.style.display = state.spectrum ? "none" : "grid";
}

function render() {
  updateModeButtons();
  renderAtlasControls();
  renderSummary();
  renderSelectionLabel();
  renderRegionsTable();
  renderStatus();
  renderCursorReadout();
  renderEmptyState();
  renderPlot();
}

function generateRegionId() {
  const id = `r_${String(state.regionCounter).padStart(4, "0")}`;
  state.regionCounter += 1;
  return id;
}

function replaceRegion(updatedRegion) {
  state.regions = state.regions.map((region) => (region.id === updatedRegion.id ? updatedRegion : region));
}

function removeRegion(regionId) {
  state.regions = state.regions.filter((region) => region.id !== regionId);
  if (state.selectedRegionId === regionId) {
    state.selectedRegionId = null;
  }
}

function readFile(file) {
  if (file && file.name.toLowerCase().endsWith(".json")) {
    readRegionsFile(file);
    return;
  }
  readSpectrumFile(file);
}

function readSpectrumFile(file) {
  if (!file) {
    return;
  }
  const reader = new FileReader();
  reader.onload = () => {
    const { spectrum, warnings } = parseSpectrumText(String(reader.result), file.name);
    if (!spectrum) {
      setStatus(`Failed to parse ${file.name}.`, warnings);
      render();
      return;
    }

    state.spectrum = spectrum;
    state.regionSourceName = spectrum.sourceName;
    state.hoverWavelength = null;
    resetView();
    setStatus(`Loaded ${file.name} with ${spectrum.pointCount} valid points.`, warnings);
    render();
  };
  reader.onerror = () => {
    setStatus(`Could not read ${file.name}.`);
    render();
  };
  reader.readAsText(file);
}

function readRegionsFile(file) {
  if (!file) {
    return;
  }

  const reader = new FileReader();
  reader.onload = () => {
    let imported;
    try {
      imported = parseRegionsJson(String(reader.result));
    } catch (error) {
      const message = error instanceof Error ? error.message : "Could not parse region JSON.";
      setStatus(`Failed to parse ${file.name}.`, [message]);
      render();
      return;
    }

    const importedRegions = imported.regions.map((region) =>
      createRegion(generateRegionId(), region.type, region.start, region.end)
    );
    state.regions = sortRegions(importedRegions);
    state.selectedRegionId = state.regions.length > 0 ? state.regions[0].id : null;
    if (imported.sourceName) {
      state.regionSourceName = imported.sourceName;
    }

    const warnings = [];
    if (!state.spectrum) {
      warnings.push("Load a spectrum to visualize the imported regions.");
    }
    setStatus(
      `Loaded ${state.regions.length} region${state.regions.length === 1 ? "" : "s"} from ${file.name}.`,
      warnings
    );
    render();
  };
  reader.onerror = () => {
    setStatus(`Could not read ${file.name}.`);
    render();
  };
  reader.readAsText(file);
}

function getRelativeCanvasPoint(event) {
  const rect = elements.plotCanvas.getBoundingClientRect();
  return {
    x: event.clientX - rect.left,
    y: event.clientY - rect.top
  };
}

function clampDataX(value) {
  if (!state.spectrum) {
    return value;
  }
  return Math.min(state.spectrum.wavelengthMax, Math.max(state.spectrum.wavelengthMin, value));
}

function pointerToDataX(event) {
  const point = getRelativeCanvasPoint(event);
  return clampDataX(plot.screenXToDataX(point.x, state.view));
}

function beginCreateRegion(dataX) {
  state.selectedRegionId = null;
  state.interaction = {
    mode: "creating",
    regionId: null,
    anchorWavelength: dataX,
    startWavelength: dataX,
    originalStart: null,
    originalEnd: null
  };
  state.draftRegion = createRegion("draft", state.activeToolType, dataX, dataX);
}

function beginZoomBox(dataX) {
  state.interaction = {
    mode: "zooming-box",
    regionId: null,
    anchorWavelength: dataX,
    startWavelength: dataX,
    originalStart: null,
    originalEnd: null
  };
  state.zoomBox = normalizeBounds(dataX, dataX);
}

function applyZoomSpan(start, end) {
  if (!state.spectrum || !state.view) {
    return;
  }

  const bounds = normalizeBounds(clampDataX(start), clampDataX(end));
  const minWidth = getMinRegionWidthData();

  if (bounds.end - bounds.start < minWidth) {
    return;
  }

  updateViewForXRange(bounds.start, bounds.end);
}

function zoomByFactor(factor) {
  if (!state.spectrum || !state.view) {
    return;
  }

  const currentSpan = state.view.xMax - state.view.xMin;
  const center = state.view.xMin + currentSpan / 2;
  const fullSpan = state.spectrum.wavelengthMax - state.spectrum.wavelengthMin;
  let nextSpan = currentSpan * factor;

  nextSpan = Math.max(fullSpan * 0.001, Math.min(fullSpan, nextSpan));

  let xMin = center - nextSpan / 2;
  let xMax = center + nextSpan / 2;

  if (xMin < state.spectrum.wavelengthMin) {
    xMin = state.spectrum.wavelengthMin;
    xMax = xMin + nextSpan;
  }
  if (xMax > state.spectrum.wavelengthMax) {
    xMax = state.spectrum.wavelengthMax;
    xMin = xMax - nextSpan;
  }

  updateViewForXRange(xMin, xMax);
}

function zoomToSelectedRegion() {
  const region = getSelectedRegion();
  if (!region || !state.spectrum) {
    setStatus("No region selected for zoom.", state.warnings);
    render();
    return;
  }

  const width = region.end - region.start;
  const padding = Math.max(width * 0.08, getMinRegionWidthData() * 2);
  applyZoomSpan(region.start - padding, region.end + padding);
  setStatus(`Zoomed to selected region ${region.id}.`, state.warnings);
  render();
}

function onPointerDown(event) {
  if (!state.spectrum || !state.view) {
    return;
  }

  const point = getRelativeCanvasPoint(event);
  const plotBounds = plot.getPlotBounds();
  if (
    point.x < plotBounds.left ||
    point.x > plotBounds.right ||
    point.y < plotBounds.top ||
    point.y > plotBounds.bottom
  ) {
    return;
  }

  elements.plotCanvas.setPointerCapture(event.pointerId);

  const dataX = pointerToDataX(event);
  if (event.shiftKey) {
    state.interaction = {
      mode: "panning",
      regionId: null,
      anchorWavelength: dataX,
      startWavelength: null,
      originalStart: state.view.xMin,
      originalEnd: state.view.xMax
    };
    render();
    return;
  }

  if (state.zoomMode || event.altKey) {
    beginZoomBox(dataX);
    render();
    return;
  }

  const hit = hitTestRegions(
    state.regions,
    (value) => plot.dataXToScreenX(value, state.view),
    plotBounds,
    point.x,
    point.y,
    state.activeToolType
  );
  if (hit) {
    const region = state.regions.find((item) => item.id === hit.regionId);
    state.selectedRegionId = hit.regionId;
    state.interaction = {
      mode: hit.mode,
      regionId: hit.regionId,
      anchorWavelength: dataX,
      startWavelength: null,
      originalStart: region.start,
      originalEnd: region.end
    };
    render();
    return;
  }

  beginCreateRegion(dataX);
  render();
}

function onPointerMove(event) {
  if (!state.spectrum || !state.view) {
    return;
  }

  const point = getRelativeCanvasPoint(event);
  const plotBounds = plot.getPlotBounds();
  if (point.x >= plotBounds.left && point.x <= plotBounds.right) {
    state.hoverWavelength = clampDataX(plot.screenXToDataX(point.x, state.view));
  } else {
    state.hoverWavelength = null;
  }

  const minWidth = getMinRegionWidthData();
  const dataX = pointerToDataX(event);

  if (state.interaction.mode === "creating" && state.draftRegion) {
    const bounds = normalizeBounds(state.interaction.anchorWavelength, dataX);
    state.draftRegion = {
      ...state.draftRegion,
      start: bounds.start,
      end: bounds.end
    };
  } else if (state.interaction.mode === "zooming-box" && state.zoomBox) {
    state.zoomBox = normalizeBounds(state.interaction.anchorWavelength, dataX);
  } else if (state.interaction.mode === "moving" && state.interaction.regionId) {
    const delta = dataX - state.interaction.anchorWavelength;
    const original = {
      id: state.interaction.regionId,
      type: getSelectedRegion().type,
      start: state.interaction.originalStart + delta,
      end: state.interaction.originalEnd + delta
    };
    replaceRegion(clampRegion(original, state.spectrum.wavelengthMin, state.spectrum.wavelengthMax, minWidth));
  } else if (
    (state.interaction.mode === "resizing-left" || state.interaction.mode === "resizing-right") &&
    state.interaction.regionId
  ) {
    const region = getSelectedRegion();
    if (region) {
      replaceRegion(
        resizeRegion(
          region,
          state.interaction.mode === "resizing-left" ? "left" : "right",
          dataX,
          minWidth
        )
      );
    }
  } else if (state.interaction.mode === "panning") {
    const delta = dataX - state.interaction.anchorWavelength;
    const span = state.interaction.originalEnd - state.interaction.originalStart;
    let xMin = state.interaction.originalStart - delta;
    let xMax = xMin + span;
    const min = state.spectrum.wavelengthMin;
    const max = state.spectrum.wavelengthMax;

    if (xMin < min) {
      xMin = min;
      xMax = min + span;
    }
    if (xMax > max) {
      xMax = max;
      xMin = max - span;
    }
    updateViewForXRange(xMin, xMax);
  }

  render();
}

function onPointerUp(event) {
  if (!state.spectrum || !state.view) {
    return;
  }
  if (elements.plotCanvas.hasPointerCapture(event.pointerId)) {
    elements.plotCanvas.releasePointerCapture(event.pointerId);
  }

  if (state.interaction.mode === "creating" && state.draftRegion) {
    const minWidth = getMinRegionWidthData();
    if (regionWidth(state.draftRegion) >= minWidth) {
      const region = createRegion(
        generateRegionId(),
        state.draftRegion.type,
        state.draftRegion.start,
        state.draftRegion.end
      );
      state.regions = [...state.regions, region];
      state.selectedRegionId = region.id;
      setStatus(`Created ${region.type} region ${formatNumber(region.start)} to ${formatNumber(region.end)}.`, state.warnings);
    }
  } else if (state.interaction.mode === "zooming-box" && state.zoomBox) {
    applyZoomSpan(state.zoomBox.start, state.zoomBox.end);
    setStatus(
      `Zoomed to ${formatNumber(state.view.xMin)} to ${formatNumber(state.view.xMax)}.`,
      state.warnings
    );
  }

  state.interaction = {
    mode: null,
    regionId: null,
    anchorWavelength: null,
    startWavelength: null,
    originalStart: null,
    originalEnd: null
  };
  state.draftRegion = null;
  state.zoomBox = null;
  render();
}

function onPointerLeave() {
  state.hoverWavelength = null;
  render();
}

function onWheel(event) {
  if (!state.spectrum || !state.view) {
    return;
  }

  event.preventDefault();
  const point = getRelativeCanvasPoint(event);
  const plotBounds = plot.getPlotBounds();
  if (point.x < plotBounds.left || point.x > plotBounds.right) {
    return;
  }

  const cursorX = clampDataX(plot.screenXToDataX(point.x, state.view));
  const currentSpan = state.view.xMax - state.view.xMin;
  const zoomFactor = event.deltaY < 0 ? 0.85 : 1.18;
  let nextSpan = currentSpan * zoomFactor;
  const fullSpan = state.spectrum.wavelengthMax - state.spectrum.wavelengthMin;

  nextSpan = Math.max(fullSpan * 0.001, Math.min(fullSpan, nextSpan));
  const ratio = (cursorX - state.view.xMin) / currentSpan;

  let xMin = cursorX - ratio * nextSpan;
  let xMax = xMin + nextSpan;

  if (xMin < state.spectrum.wavelengthMin) {
    xMin = state.spectrum.wavelengthMin;
    xMax = xMin + nextSpan;
  }
  if (xMax > state.spectrum.wavelengthMax) {
    xMax = state.spectrum.wavelengthMax;
    xMin = xMax - nextSpan;
  }

  updateViewForXRange(xMin, xMax);
  render();
}

function updateRegionBoundary(regionId, field, rawValue) {
  const parsed = Number(rawValue);
  if (!Number.isFinite(parsed)) {
    setStatus("Boundary edit ignored because the value is not numeric.", state.warnings);
    render();
    return;
  }

  const region = state.regions.find((item) => item.id === regionId);
  if (!region) {
    return;
  }

  const updated = {
    ...region,
    [field]: displayToRawWavelength(parsed)
  };
  const bounds = normalizeBounds(updated.start, updated.end);
  updated.start = bounds.start;
  updated.end = bounds.end;

  const minWidth = getMinRegionWidthData() || Number.EPSILON;
  const clamped = clampRegion(updated, state.spectrum.wavelengthMin, state.spectrum.wavelengthMax, minWidth);
  replaceRegion(clamped);
  state.selectedRegionId = regionId;
  setStatus(`Updated region ${regionId}.`, state.warnings);
  render();
}

function onRegionsTableClick(event) {
  const row = event.target.closest("tr[data-region-id]");
  if (!row) {
    return;
  }

  const regionId = row.dataset.regionId;
  if (event.target.classList.contains("row-delete")) {
    removeRegion(regionId);
    setStatus(`Deleted region ${regionId}.`, state.warnings);
    render();
    return;
  }

  if (event.target.classList.contains("row-select") || row) {
    state.selectedRegionId = regionId;
    render();
  }
}

function onRegionsTableChange(event) {
  const input = event.target.closest(".boundary-input");
  if (!input) {
    return;
  }

  const row = event.target.closest("tr[data-region-id]");
  if (!row) {
    return;
  }

  updateRegionBoundary(row.dataset.regionId, input.dataset.field, input.value);
}

function deleteSelectedRegion() {
  if (!state.selectedRegionId) {
    setStatus("No region selected.", state.warnings);
    render();
    return;
  }
  const regionId = state.selectedRegionId;
  removeRegion(regionId);
  setStatus(`Deleted region ${regionId}.`, state.warnings);
  render();
}

function bindEvents() {
  ["dragenter", "dragover", "dragleave", "drop"].forEach((type) => {
    document.addEventListener(type, (event) => {
      event.preventDefault();
    });
  });

  elements.fileInput.addEventListener("change", (event) => {
    const file = event.target.files && event.target.files[0];
    readFile(file);
    event.target.value = "";
  });
  elements.regionsInput.addEventListener("change", (event) => {
    const file = event.target.files && event.target.files[0];
    readRegionsFile(file);
    event.target.value = "";
  });

  ["dragenter", "dragover"].forEach((type) => {
    elements.dropZone.addEventListener(type, (event) => {
      event.preventDefault();
      elements.dropZone.classList.add("is-over");
    });
  });

  ["dragleave", "drop"].forEach((type) => {
    elements.dropZone.addEventListener(type, (event) => {
      event.preventDefault();
      elements.dropZone.classList.remove("is-over");
    });
  });

  elements.dropZone.addEventListener("drop", (event) => {
    const file = event.dataTransfer?.files?.[0];
    readFile(file);
  });

  elements.modeInclusion.addEventListener("click", () => {
    state.activeToolType = "inclusion";
    render();
  });
  elements.modeExclusion.addEventListener("click", () => {
    state.activeToolType = "exclusion";
    render();
  });
  elements.unitAngstrom.addEventListener("click", () => {
    state.wavelengthUnit = "angstrom";
    render();
  });
  elements.unitNm.addEventListener("click", () => {
    state.wavelengthUnit = "nm";
    render();
  });
  elements.unitRaw.addEventListener("click", () => {
    state.wavelengthUnit = "raw";
    render();
  });
  elements.yScaleAuto.addEventListener("click", () => {
    state.yScaleMode = "auto";
    if (state.view) {
      updateViewForXRange(state.view.xMin, state.view.xMax);
    }
    setStatus(`Y scale set to ${Y_SCALE_LABELS[state.yScaleMode]}.`, state.warnings);
    render();
  });
  elements.yScaleZero.addEventListener("click", () => {
    state.yScaleMode = "zero";
    if (state.view) {
      updateViewForXRange(state.view.xMin, state.view.xMax);
    }
    setStatus(`Y scale set to ${Y_SCALE_LABELS[state.yScaleMode]}.`, state.warnings);
    render();
  });
  elements.yScaleFull.addEventListener("click", () => {
    state.yScaleMode = "full";
    if (state.view) {
      updateViewForXRange(state.view.xMin, state.view.xMax);
    }
    setStatus(`Y scale set to ${Y_SCALE_LABELS[state.yScaleMode]}.`, state.warnings);
    render();
  });

  elements.deleteSelected.addEventListener("click", deleteSelectedRegion);
  elements.zoomBoxToggle.addEventListener("click", () => {
    state.zoomMode = !state.zoomMode;
    render();
  });
  elements.zoomIn.addEventListener("click", () => {
    zoomByFactor(0.7);
    setStatus("Zoomed in.", state.warnings);
    render();
  });
  elements.zoomOut.addEventListener("click", () => {
    zoomByFactor(1.4);
    setStatus("Zoomed out.", state.warnings);
    render();
  });
  elements.zoomSelected.addEventListener("click", zoomToSelectedRegion);
  elements.atomicAtlasSelect.addEventListener("change", (event) => {
    const nextAtlasId = event.target.value;
    const previousAtlasId = state.selectedAtlasId;
    setActiveAtlas(nextAtlasId);
    if (state.selectedAtlasId !== previousAtlasId) {
      setStatus(getAtomicLabelsStatusMessage(), state.warnings);
    }
    render();
  });
  elements.resetView.addEventListener("click", () => {
    resetView();
    state.zoomBox = null;
    setStatus(`Reset view with ${Y_SCALE_LABELS[state.yScaleMode]} y scale.`, state.warnings);
    render();
  });

  elements.exportRegions.addEventListener("click", () => {
    const message = exportRegions(
      state.regions,
      state.spectrum?.sourceName || state.regionSourceName || "regions"
    );
    setStatus(message, state.warnings);
    render();
  });

  elements.regionsBody.addEventListener("click", onRegionsTableClick);
  elements.regionsBody.addEventListener("change", onRegionsTableChange);

  elements.plotCanvas.addEventListener("pointerdown", onPointerDown);
  elements.plotCanvas.addEventListener("pointermove", onPointerMove);
  elements.plotCanvas.addEventListener("pointerup", onPointerUp);
  elements.plotCanvas.addEventListener("pointercancel", onPointerUp);
  elements.plotCanvas.addEventListener("pointerleave", onPointerLeave);
  elements.plotCanvas.addEventListener("wheel", onWheel, { passive: false });

  window.addEventListener("keydown", (event) => {
    const target = event.target;
    const tagName = target && "tagName" in target ? String(target.tagName).toLowerCase() : "";
    if (tagName === "input" || tagName === "textarea") {
      return;
    }
    if (event.key === "Delete") {
      deleteSelectedRegion();
    }
  });

  window.addEventListener("resize", () => {
    plot.resize();
    render();
  });
}

function init() {
  elements.atomicAtlasSelect.innerHTML = [
    `<option value="${NO_ATLAS_ID}">No atomic labels</option>`,
    ...ATLAS_CATALOG.map((atlas) => `<option value="${atlas.id}">Atlas: ${atlas.label}</option>`)
  ].join("");
  setActiveAtlas(DEFAULT_ATLAS_ID);
  plot.resize();
  bindEvents();
  render();
}

init();
