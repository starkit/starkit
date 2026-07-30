const HANDLE_SIZE_PX = 7;

export function normalizeBounds(start, end) {
  if (start <= end) {
    return { start, end };
  }
  return { start: end, end: start };
}

export function createRegion(id, type, start, end) {
  const bounds = normalizeBounds(start, end);
  return {
    id,
    type,
    start: bounds.start,
    end: bounds.end
  };
}

export function sortRegions(regions) {
  return [...regions].sort((a, b) => {
    if (a.start !== b.start) {
      return a.start - b.start;
    }
    return a.end - b.end;
  });
}

export function regionWidth(region) {
  return region.end - region.start;
}

export function clampRegion(region, minX, maxX, minWidth) {
  let start = region.start;
  let end = region.end;

  if (start < minX) {
    end += minX - start;
    start = minX;
  }
  if (end > maxX) {
    start -= end - maxX;
    end = maxX;
  }

  start = Math.max(minX, start);
  end = Math.min(maxX, end);

  if (end - start < minWidth) {
    end = Math.min(maxX, start + minWidth);
    start = Math.max(minX, end - minWidth);
  }

  return { ...region, start, end };
}

export function resizeRegion(region, edge, value, minWidth) {
  let start = region.start;
  let end = region.end;

  if (edge === "left") {
    start = Math.min(value, end - minWidth);
  } else {
    end = Math.max(value, start + minWidth);
  }

  return { ...region, start, end };
}

export function hitTestRegions(regions, dataXToScreenX, plotBounds, x, y, allowedType = null) {
  if (
    x < plotBounds.left ||
    x > plotBounds.right ||
    y < plotBounds.top ||
    y > plotBounds.bottom
  ) {
    return null;
  }

  for (let index = regions.length - 1; index >= 0; index -= 1) {
    const region = regions[index];
    if (allowedType && region.type !== allowedType) {
      continue;
    }
    const x1 = dataXToScreenX(region.start);
    const x2 = dataXToScreenX(region.end);
    const left = Math.min(x1, x2);
    const right = Math.max(x1, x2);

    if (Math.abs(x - left) <= HANDLE_SIZE_PX) {
      return { regionId: region.id, mode: "resizing-left" };
    }
    if (Math.abs(x - right) <= HANDLE_SIZE_PX) {
      return { regionId: region.id, mode: "resizing-right" };
    }
    if (x >= left && x <= right) {
      return { regionId: region.id, mode: "moving" };
    }
  }

  return null;
}
