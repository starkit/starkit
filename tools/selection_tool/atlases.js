import { ATOMIC_LABELS } from "./atomic_labels.js";
import { ARCTURUS_ATOMIC_LABELS } from "./arcturus_labels.js";

function createAtlasDefinition(id, label, labels) {
  const sortedLabels = [...labels].sort((a, b) => a.x - b.x);
  return {
    id,
    label,
    labels: sortedLabels,
    wavelengthMin: sortedLabels[0]?.x ?? null,
    wavelengthMax: sortedLabels[sortedLabels.length - 1]?.x ?? null
  };
}

export const ATLAS_CATALOG = [
  createAtlasDefinition("cool_giant", "Cool Giant", ATOMIC_LABELS),
  createAtlasDefinition("arcturus", "Arcturus", ARCTURUS_ATOMIC_LABELS)
];

export const DEFAULT_ATLAS_ID = ATLAS_CATALOG[0].id;

export function getAtlasById(id) {
  return ATLAS_CATALOG.find((atlas) => atlas.id === id) || ATLAS_CATALOG[0];
}
