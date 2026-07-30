function isFiniteNumber(value) {
  return Number.isFinite(value);
}

export function parseSpectrumText(text, sourceName = "spectrum.txt") {
  const lines = text.split(/\r?\n/);
  const points = [];
  let skippedBlankOrComment = 0;
  let skippedInvalid = 0;
  let hadUnsortedInput = false;

  for (const rawLine of lines) {
    const line = rawLine.trim();
    if (!line || line.startsWith("#")) {
      skippedBlankOrComment += 1;
      continue;
    }

    const fields = line.split(/[\s,]+/).filter(Boolean);
    if (fields.length < 2) {
      skippedInvalid += 1;
      continue;
    }

    const wavelength = Number(fields[0]);
    const flux = Number(fields[1]);
    if (!isFiniteNumber(wavelength) || !isFiniteNumber(flux)) {
      skippedInvalid += 1;
      continue;
    }

    const lastPoint = points[points.length - 1];
    if (lastPoint && wavelength < lastPoint.wavelength) {
      hadUnsortedInput = true;
    }

    points.push({ wavelength, flux });
  }

  if (points.length === 0) {
    return {
      spectrum: null,
      warnings: [`No valid rows found in ${sourceName}.`],
      stats: { skippedBlankOrComment, skippedInvalid, validRows: 0 }
    };
  }

  const sortedPoints = [...points].sort((a, b) => a.wavelength - b.wavelength);
  const wavelengthMin = sortedPoints[0].wavelength;
  const wavelengthMax = sortedPoints[sortedPoints.length - 1].wavelength;

  let fluxMin = sortedPoints[0].flux;
  let fluxMax = sortedPoints[0].flux;
  for (const point of sortedPoints) {
    if (point.flux < fluxMin) {
      fluxMin = point.flux;
    }
    if (point.flux > fluxMax) {
      fluxMax = point.flux;
    }
  }

  if (fluxMin === fluxMax) {
    fluxMin -= 1;
    fluxMax += 1;
  }

  const warnings = [];
  if (hadUnsortedInput) {
    warnings.push("Input wavelengths were not sorted; rows were reordered by wavelength.");
  }
  if (skippedInvalid > 0) {
    warnings.push(`Skipped ${skippedInvalid} invalid data line${skippedInvalid === 1 ? "" : "s"}.`);
  }

  return {
    spectrum: {
      sourceName,
      points: sortedPoints,
      wavelengthMin,
      wavelengthMax,
      fluxMin,
      fluxMax,
      pointCount: sortedPoints.length
    },
    warnings,
    stats: {
      skippedBlankOrComment,
      skippedInvalid,
      validRows: sortedPoints.length
    }
  };
}
