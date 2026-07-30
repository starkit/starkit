const MARGINS = { top: 72, right: 18, bottom: 42, left: 72 };
const MIN_LABEL_SPACING_PX = 18;

function formatTick(value) {
  if (!Number.isFinite(value)) {
    return "";
  }
  const abs = Math.abs(value);
  if (abs >= 1000 || (abs > 0 && abs < 0.01)) {
    return value.toExponential(2);
  }
  if (abs >= 100) {
    return value.toFixed(2);
  }
  if (abs >= 1) {
    return value.toFixed(4).replace(/\.?0+$/, "");
  }
  return value.toPrecision(4);
}

function drawLine(ctx, points) {
  if (points.length === 0) {
    return;
  }

  ctx.beginPath();
  ctx.moveTo(points[0].x, points[0].y);
  for (let index = 1; index < points.length; index += 1) {
    ctx.lineTo(points[index].x, points[index].y);
  }
  ctx.stroke();
}

function findFirstLabelAtOrAfter(labels, wavelength) {
  let low = 0;
  let high = labels.length;

  while (low < high) {
    const mid = low + Math.floor((high - low) / 2);
    if (labels[mid].x < wavelength) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  return low;
}

export class SpectrumPlot {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.width = 0;
    this.height = 0;
    this.pixelRatio = window.devicePixelRatio || 1;
    this.plotBounds = { left: 0, right: 0, top: 0, bottom: 0 };
  }

  resize() {
    const rect = this.canvas.getBoundingClientRect();
    const width = Math.max(300, Math.floor(rect.width));
    const height = Math.max(260, Math.floor(rect.height));
    this.pixelRatio = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(width * this.pixelRatio);
    this.canvas.height = Math.round(height * this.pixelRatio);
    this.width = width;
    this.height = height;
    this.ctx.setTransform(this.pixelRatio, 0, 0, this.pixelRatio, 0, 0);
    this.plotBounds = {
      left: MARGINS.left,
      right: width - MARGINS.right,
      top: MARGINS.top,
      bottom: height - MARGINS.bottom
    };
  }

  getPlotBounds() {
    return this.plotBounds;
  }

  dataXToScreenX(x, view) {
    const { left, right } = this.plotBounds;
    return left + ((x - view.xMin) / (view.xMax - view.xMin)) * (right - left);
  }

  dataYToScreenY(y, view) {
    const { top, bottom } = this.plotBounds;
    return bottom - ((y - view.yMin) / (view.yMax - view.yMin)) * (bottom - top);
  }

  screenXToDataX(x, view) {
    const { left, right } = this.plotBounds;
    const ratio = (x - left) / (right - left);
    return view.xMin + ratio * (view.xMax - view.xMin);
  }

  draw({
    spectrum,
    regions,
    selectedRegionId,
    view,
    draftRegion,
    zoomBox,
    hoverWavelength,
    atomicLabels,
    showAtomicLabels,
    wavelengthFormatter,
    wavelengthUnitLabel
  }) {
    const ctx = this.ctx;
    ctx.clearRect(0, 0, this.width, this.height);

    this.drawBackground();
    if (!spectrum || !view) {
      return;
    }

    this.drawAxes(view, wavelengthFormatter, wavelengthUnitLabel);
    this.drawRegions(regions, view, selectedRegionId);
    if (draftRegion) {
      this.drawRegions([draftRegion], view, draftRegion.id);
    }
    if (zoomBox) {
      this.drawZoomBox(zoomBox, view);
    }
    this.drawSpectrum(spectrum, view);
    if (showAtomicLabels && atomicLabels && atomicLabels.length > 0) {
      this.drawAtomicLabels(atomicLabels, view);
    }
    if (hoverWavelength !== null) {
      this.drawCursor(hoverWavelength, view);
    }
  }

  drawBackground() {
    const ctx = this.ctx;
    const { left, top, right, bottom } = this.plotBounds;
    ctx.fillStyle = "#fffef9";
    ctx.fillRect(0, 0, this.width, this.height);
    ctx.fillStyle = "#fffdfa";
    ctx.fillRect(left, top, right - left, bottom - top);
  }

  drawAxes(view, wavelengthFormatter, wavelengthUnitLabel) {
    const ctx = this.ctx;
    const { left, top, right, bottom } = this.plotBounds;
    const xTicks = 6;
    const yTicks = 5;

    ctx.strokeStyle = "#dcd5c6";
    ctx.fillStyle = "#6d665a";
    ctx.font = '12px "Segoe UI", sans-serif';
    ctx.lineWidth = 1;

    for (let tick = 0; tick <= xTicks; tick += 1) {
      const ratio = tick / xTicks;
      const x = left + ratio * (right - left);
      const value = view.xMin + ratio * (view.xMax - view.xMin);

      ctx.beginPath();
      ctx.moveTo(x, top);
      ctx.lineTo(x, bottom);
      ctx.stroke();

      const label = wavelengthFormatter ? wavelengthFormatter(value) : formatTick(value);
      ctx.fillText(label, x - 24, bottom + 20);
    }

    for (let tick = 0; tick <= yTicks; tick += 1) {
      const ratio = tick / yTicks;
      const y = bottom - ratio * (bottom - top);
      const value = view.yMin + ratio * (view.yMax - view.yMin);

      ctx.beginPath();
      ctx.moveTo(left, y);
      ctx.lineTo(right, y);
      ctx.stroke();

      ctx.fillText(formatTick(value), 10, y + 4);
    }

    ctx.strokeStyle = "#847b6b";
    ctx.beginPath();
    ctx.moveTo(left, top);
    ctx.lineTo(left, bottom);
    ctx.lineTo(right, bottom);
    ctx.stroke();

    ctx.fillStyle = "#6d665a";
    ctx.font = '12px "Segoe UI", sans-serif';
    ctx.fillText(`Wavelength (${wavelengthUnitLabel || "raw"})`, right - 142, this.height - 8);
    ctx.save();
    ctx.translate(16, top + (bottom - top) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("Flux", 0, 0);
    ctx.restore();
  }

  drawSpectrum(spectrum, view) {
    const ctx = this.ctx;
    ctx.save();
    ctx.beginPath();
    ctx.rect(
      this.plotBounds.left,
      this.plotBounds.top,
      this.plotBounds.right - this.plotBounds.left,
      this.plotBounds.bottom - this.plotBounds.top
    );
    ctx.clip();

    ctx.strokeStyle = "#214f7a";
    ctx.lineWidth = 1.1;
    drawLine(
      ctx,
      spectrum.points.map((point) => ({
        x: this.dataXToScreenX(point.wavelength, view),
        y: this.dataYToScreenY(point.flux, view)
      }))
    );
    ctx.restore();
  }

  drawRegions(regions, view, selectedRegionId) {
    const ctx = this.ctx;
    const { top, bottom } = this.plotBounds;

    ctx.save();
    ctx.beginPath();
    ctx.rect(
      this.plotBounds.left,
      this.plotBounds.top,
      this.plotBounds.right - this.plotBounds.left,
      this.plotBounds.bottom - this.plotBounds.top
    );
    ctx.clip();

    for (const region of regions) {
      const left = this.dataXToScreenX(region.start, view);
      const right = this.dataXToScreenX(region.end, view);
      const stroke = region.type === "inclusion" ? "#14553d" : "#8f2c1d";
      const fill = region.type === "inclusion" ? "rgba(31, 111, 80, 0.2)" : "rgba(165, 48, 30, 0.2)";

      ctx.fillStyle = fill;
      ctx.strokeStyle = region.id === selectedRegionId ? "#d39d1a" : stroke;
      ctx.lineWidth = region.id === selectedRegionId ? 2.5 : 1.5;

      ctx.fillRect(left, top, right - left, bottom - top);
      ctx.strokeRect(left, top, right - left, bottom - top);

      ctx.fillStyle = ctx.strokeStyle;
      ctx.fillRect(left - 1, top, 3, bottom - top);
      ctx.fillRect(right - 1, top, 3, bottom - top);
    }

    ctx.restore();
  }

  drawCursor(wavelength, view) {
    const ctx = this.ctx;
    const { top, bottom } = this.plotBounds;
    const x = this.dataXToScreenX(wavelength, view);
    if (x < this.plotBounds.left || x > this.plotBounds.right) {
      return;
    }

    ctx.save();
    ctx.strokeStyle = "rgba(33, 79, 122, 0.55)";
    ctx.setLineDash([5, 5]);
    ctx.beginPath();
    ctx.moveTo(x, top);
    ctx.lineTo(x, bottom);
    ctx.stroke();
    ctx.restore();
  }

  drawZoomBox(zoomBox, view) {
    const ctx = this.ctx;
    const { top, bottom } = this.plotBounds;
    const left = this.dataXToScreenX(zoomBox.start, view);
    const right = this.dataXToScreenX(zoomBox.end, view);

    ctx.save();
    ctx.fillStyle = "rgba(33, 79, 122, 0.14)";
    ctx.strokeStyle = "#214f7a";
    ctx.lineWidth = 1.5;
    ctx.setLineDash([6, 4]);
    ctx.fillRect(left, top, right - left, bottom - top);
    ctx.strokeRect(left, top, right - left, bottom - top);
    ctx.restore();
  }

  drawAtomicLabels(labels, view) {
    const ctx = this.ctx;
    const topGuideY = this.plotBounds.top - 8;
    const labelBaseY = this.plotBounds.top - 12;
    const visible = [];
    let lastAcceptedX = -Infinity;
    const startIndex = findFirstLabelAtOrAfter(labels, view.xMin);

    for (let index = startIndex; index < labels.length; index += 1) {
      const label = labels[index];
      if (label.x > view.xMax) {
        break;
      }
      const screenX = this.dataXToScreenX(label.x, view);
      if (screenX - lastAcceptedX < MIN_LABEL_SPACING_PX) {
        continue;
      }
      visible.push({ ...label, screenX });
      lastAcceptedX = screenX;
    }

    ctx.save();
    ctx.strokeStyle = "rgba(89, 76, 45, 0.45)";
    ctx.fillStyle = "#594c2d";
    ctx.font = '11px "Segoe UI", sans-serif';
    ctx.lineWidth = 1;

    for (const label of visible) {
      ctx.beginPath();
      ctx.moveTo(label.screenX, this.plotBounds.top);
      ctx.lineTo(label.screenX, topGuideY);
      ctx.stroke();

      ctx.save();
      ctx.translate(label.screenX + 3, labelBaseY);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(label.text, 0, 0);
      ctx.restore();
    }

    ctx.restore();
  }
}
