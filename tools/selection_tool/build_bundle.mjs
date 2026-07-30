import { readFileSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const directory = dirname(fileURLToPath(import.meta.url));
const inputFiles = [
  "parser.js",
  "export.js",
  "atomic_labels.js",
  "arcturus_labels.js",
  "atlases.js",
  "regions.js",
  "plot.js",
  "app.js"
];
const standaloneOutputFile = join(directory, "selection_tool.html");
const htmlTemplateFile = join(directory, "index.html");
const stylesFile = join(directory, "styles.css");

function transformModule(source, fileName) {
  const transformed = source
    .replace(/^import\s+[^;]+;\n/gm, "")
    .replace(/^export\s+class\s+/gm, "class ")
    .replace(/^export\s+function\s+/gm, "function ")
    .replace(/^export\s+const\s+/gm, "const ");

  if (/^\s*export\s+/m.test(transformed)) {
    throw new Error(`Unsupported export syntax remains in ${fileName}.`);
  }

  return `// ${fileName}\n${transformed.trim()}\n`;
}

const bundleParts = ['"use strict";\n'];

for (const fileName of inputFiles) {
  const source = readFileSync(join(directory, fileName), "utf8");
  bundleParts.push(transformModule(source, fileName));
}

const bundleSource = bundleParts.join("\n");
const htmlTemplate = readFileSync(htmlTemplateFile, "utf8");
const stylesSource = readFileSync(stylesFile, "utf8");
const escapedBundleSource = bundleSource.replace(/<\/script/gi, "<\\/script");

const standaloneHtml = htmlTemplate
  .replace('<link rel="stylesheet" href="./styles.css" />', `<style>\n${stylesSource}\n    </style>`)
  .replace('<script type="module" src="./app.js"></script>', `<script>\n${escapedBundleSource}\n    </script>`);

writeFileSync(standaloneOutputFile, standaloneHtml, "utf8");
