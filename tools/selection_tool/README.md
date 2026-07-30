# Selection Tool

This folder contains a small standalone browser tool for loading a spectrum, defining inclusion/exclusion regions, and exporting those regions as JSON.

## Important note about `file:///`

Opening the tool directly from disk can break ES module loading in Chromium-based browsers. For that reason, the distributable version is a generated file that does not rely on module imports at runtime.

## Source files

Edit these files when changing the tool:

- `app.js`
- `atlases.js`
- `plot.js`
- `regions.js`
- `parser.js`
- `export.js`
- `atomic_labels.js`
- `arcturus_labels.js`
- `generate_arcturus_labels.mjs`
- `styles.css`
- `index.html`

## Generated files

This file is the build output:

- `selection_tool.html`
  Single-file distribution with CSS and JS inlined. This is the file to send to someone else and the file to open directly via `file:///`.

## Rebuild after edits

If any source file changes, regenerate the standalone distribution with:

```bash
node tools/selection_tool/build_bundle.mjs
```

That command rebuilds `selection_tool.html`.

If the Arcturus source file changes, regenerate `arcturus_labels.js` first:

```bash
node tools/selection_tool/generate_arcturus_labels.mjs \
  /path/to/arcturus-atlas_line_ids_visual_ir_atoms.txt
```

## Distribution options

Send `selection_tool.html`.

## Practical workflow

1. Edit the source files.
2. Run `node tools/selection_tool/build_bundle.mjs`.
3. Open `selection_tool.html` or refresh it in the browser.

Closing and reopening the browser does not require a rebuild by itself. Rebuild only when the source files have changed.

## About `index.html`

`index.html` is the source template for the generated standalone file. It loads `app.js` as an ES module and may require a local web server if opened directly in some browsers. For local file usage, use `selection_tool.html` instead.
