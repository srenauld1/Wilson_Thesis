import fs from 'node:fs/promises';
import {FileBlob,SpreadsheetFile} from '@oai/artifact-tool';
const w=await SpreadsheetFile.importXlsx(await FileBlob.load('/Users/sophiarenauld/stacks/pizza_talk_analysis/good_behavior_fly_inventory.xlsx'));
console.log((await w.inspect({kind:'table',range:"'Trial records'!A45:J53",tableMaxRows:9,tableMaxCols:10,maxChars:2500})).ndjson);
const b=await w.render({sheetName:'Trial records',range:'A45:F53',scale:1.5});await fs.writeFile('outputs/pizza_heat_20261028/before.png',new Uint8Array(await b.arrayBuffer()));
