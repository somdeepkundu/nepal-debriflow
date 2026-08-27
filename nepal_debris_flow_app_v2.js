/*******************************************************************************
 * NEPAL DEBRIS FLOW — DOWNSTREAM EXPOSURE, IMAGERY & RAINFALL ANOMALY
 * Rasuwa River Valley (28.279 N, 85.378 E) — event 2026-08-25
 *
 * Imagery   : pre-event 1-month median composite vs post-event scene (Aug 27)
 * Hydrology : downstream-only trace via RiverATLAS NEXT_DOWN, clipped to AOI
 * Rainfall  : CHIRPS cumulative, event year against the prior-year normal
 *
 * Theme: Ocean Depths
 * Somdeep Kundu · RuDRA Lab · C-TARA, IIT Bombay
 ******************************************************************************/

// ============================================================================
// 0. THEME — Ocean Depths
// ============================================================================
var T = {
  navy:    '#1a2332',
  navy2:   '#243447',
  teal:    '#2d8b8b',
  seafoam: '#a8dadc',
  cream:   '#f1faee',
  alert:   '#e63946',
  border:  '#c9d6d6',
  textSub: '#5c6b73'
};

// ============================================================================
// 1. CONFIG
// ============================================================================
var CFG = {
  lon: 85.378,
  lat: 28.279,
  zoom: 13,

  eventDate: '2026-08-25',
  postAssetId: 'COPERNICUS/S2_SR_HARMONIZED/' +
               '20260827T045659_20260827T051017_T45RUM',
  baselineDays: 30,     // pre-event composite window
  csThreshold: 0.6,     // Cloud Score+ cutoff

  basinLevel: 10,
  aoiRadiusM: 10000,    // study area — analysis is confined to this
  corridorM:  100,
  snapM:      600,
  reachHops:  80,
  searchM:    15000,
  demGuardM:  20,

  iso3: 'NPL',
  chirpsScale: 5566,

  years: [2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026],
  eventYear: 2026,
  monthsAvailable: 8
};

var EVENT = ee.Date(CFG.eventDate);
var pt = ee.Geometry.Point([CFG.lon, CFG.lat]);
var aoi = pt.buffer(CFG.aoiRadiusM);
var searchArea = pt.buffer(CFG.searchM);

// ============================================================================
// 2. IMAGERY — pre-event composite & post-event scene
// ============================================================================
var preCollection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(EVENT.advance(-CFG.baselineDays, 'day'), EVENT)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80));

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var preMasked = preCollection.linkCollection(csPlus, ['cs'])
  .map(function (img) {
    return img.updateMask(img.select('cs').gte(CFG.csThreshold));
  });

var preImage  = preMasked.median().clip(aoi);
var postImage = ee.Image(CFG.postAssetId).clip(aoi);

// Band-combination profiles
var vizModes = {
  'SWIR FCC [B11, B8, B4]': {
    bands: ['B11', 'B8', 'B4'], min: 150, max: 3500, gamma: 1.2, key: 'swir'
  },
  'Standard FCC [B8, B4, B3]': {
    bands: ['B8', 'B4', 'B3'], min: 200, max: 4500, gamma: 1.2, key: 'fcc'
  },
  'Natural Colour [B4, B3, B2]': {
    bands: ['B4', 'B3', 'B2'], min: 150, max: 3000, gamma: 1.3, key: 'tcc'
  }
};
var defaultMode = 'SWIR FCC [B11, B8, B4]';

// Persistent layer handles so the dropdown can restyle without a rebuild
var layerPre  = ui.Map.Layer(preImage,  vizModes[defaultMode],
                             'Pre-event composite', false);
var layerPost = ui.Map.Layer(postImage, vizModes[defaultMode],
                             'Post-event (Aug 27)', true);

// ============================================================================
// 3. DOWNSTREAM RIVER TRACE
// ============================================================================
var riverAtlas = ee.FeatureCollection(
  'projects/sat-io/open-datasets/HydroAtlas/RiverAtlas_v10');
var riversRegion = riverAtlas.filterBounds(searchArea);

var seed = ee.Feature(riversRegion.filterBounds(pt.buffer(CFG.snapM))
             .map(function (f) {
               return f.set('_d', f.geometry().distance(pt, 10));
             })
             .sort('_d').first());

var seedId = ee.Number(seed.get('HYRIV_ID'));

// Walk NEXT_DOWN hop by hop. NEXT_DOWN = 0 marks a terminal reach, so the
// chain pads with zeros once it runs out; those are dropped afterwards.
var chain = ee.List(
  ee.List.sequence(1, CFG.reachHops).iterate(function (i, acc) {
    acc = ee.List(acc);
    var cur = ee.Number(acc.get(-1));
    var f = riversRegion.filter(ee.Filter.eq('HYRIV_ID', cur)).first();
    var nxt = ee.Algorithms.If(f, ee.Feature(f).get('NEXT_DOWN'), 0);
    return acc.add(ee.Number(nxt));
  }, ee.List([seedId]))
).distinct().removeAll([0]);

var downstreamRaw = riversRegion.filter(ee.Filter.inList('HYRIV_ID', chain));

// AW3D30 is an ImageCollection — it must be mosaicked before use
var dem = ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1').mosaic().select('DSM');

var seedElev = ee.Number(dem.reduceRegion({
  reducer: ee.Reducer.first(), geometry: pt, scale: 30
}).values().get(0));

// Confine to the study area, then reject anything above the epicentre
var downstream = downstreamRaw
  .filterBounds(aoi)
  .map(function (f) {
    var e = dem.reduceRegion({
      reducer: ee.Reducer.min(), geometry: f.geometry(), scale: 90,
      bestEffort: true
    }).values().get(0);
    return f.set('_elev', e);
  })
  .filter(ee.Filter.lte('_elev', seedElev.add(CFG.demGuardM)));

var upstream = riversRegion.filterBounds(aoi)
                 .filter(ee.Filter.inList('HYRIV_ID', chain).not());

var corridor = downstream.geometry(10)
                 .buffer(CFG.corridorM, 10)
                 .intersection(aoi, 10);

// ============================================================================
// 4. CATCHMENT & BUILDINGS
// ============================================================================
var basin = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_' + CFG.basinLevel)
              .filterBounds(pt);
var basinGeom = basin.geometry();

var buildings = ee.FeatureCollection(
  'projects/sat-io/open-datasets/VIDA_COMBINED/' + CFG.iso3)
  .filterBounds(aoi);

// Fallback if VIDA has no NPL tile:
// var buildings = ee.FeatureCollection('GOOGLE/Research/open-buildings/v3/polygons')
//                   .filterBounds(aoi);

var exposed = buildings.filterBounds(corridor);

// ============================================================================
// 5. RAINFALL — CHIRPS monthly + running cumulative, per year
// ============================================================================
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD');
var months = ee.List.sequence(1, 12);

var bandNames = months.map(function (m) {
    return ee.String('m').cat(ee.Number(m).format('%02d'));
  }).cat(months.map(function (m) {
    return ee.String('c').cat(ee.Number(m).format('%02d'));
  }));

// An empty date window returns a BANDLESS image, which breaks the rename
// below. CHIRPS has no data for the tail of the event year, so guard it.
function safeSum(coll) {
  return ee.Image(ee.Algorithms.If(coll.size().gt(0), coll.sum(), ee.Image(0)))
           .rename('p');
}

function yearStack(y) {
  var jan1 = ee.Date.fromYMD(y, 1, 1);
  var monthly = months.map(function (m) {
    var s = ee.Date.fromYMD(y, ee.Number(m), 1);
    return safeSum(chirps.filterDate(s, s.advance(1, 'month')));
  });
  var cumulative = months.map(function (m) {
    return safeSum(chirps.filterDate(jan1, jan1.advance(ee.Number(m), 'month')));
  });
  return ee.ImageCollection.fromImages(monthly.cat(cumulative))
           .toBands().rename(bandNames);
}

function yearStats(y) {
  return yearStack(y).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: basinGeom,
    scale: CFG.chirpsScale, maxPixels: 1e9, bestEffort: true, tileScale: 2
  });
}

var stormColl = chirps.filterDate(EVENT.advance(-45, 'day'),
                                  EVENT.advance(15, 'day'));

// ============================================================================
// 6. UI SHELL — Ocean Depths
// ============================================================================
Map.setOptions('SATELLITE');
Map.setCenter(CFG.lon, CFG.lat, CFG.zoom);
Map.style().set('cursor', 'crosshair');

var panel = ui.Panel({style: {width: '400px', padding: '0',
                              backgroundColor: T.cream}});
ui.root.insert(0, panel);

var header = ui.Panel({style: {backgroundColor: T.navy,
                               padding: '16px 16px 14px 16px',
                               margin: '0 0 10px 0'}});
header.add(ui.Label('Downstream Exposure & Rainfall', {
  fontSize: '17px', fontWeight: 'bold', color: T.cream,
  backgroundColor: T.navy, margin: '0 0 3px 0'
}));
header.add(ui.Label('Rasuwa River Valley · debris flow ' + CFG.eventDate, {
  fontSize: '11.5px', color: T.seafoam, backgroundColor: T.navy, margin: '0'
}));
panel.add(header);

var body = ui.Panel({style: {padding: '0 10px 14px 10px',
                             backgroundColor: T.cream}});
panel.add(body);

function makeKpi(label) {
  var lab = ui.Label(label, {fontSize: '10px', color: T.textSub,
                             backgroundColor: '#ffffff', margin: '7px 0 0 9px'});
  var val = ui.Label('…', {fontSize: '20px', fontWeight: 'bold', color: T.teal,
                           backgroundColor: '#ffffff', margin: '1px 0 0 9px'});
  var sub = ui.Label(' ', {fontSize: '10px', color: T.textSub,
                           backgroundColor: '#ffffff', margin: '0 0 7px 9px'});
  return {panel: ui.Panel([lab, val, sub], null, {
            width: '180px', margin: '3px', backgroundColor: '#ffffff',
            border: '1px solid ' + T.border}),
          title: lab, value: val, sub: sub};
}

var kReach   = makeKpi('DOWNSTREAM REACH');
var kBasin   = makeKpi('CATCHMENT AREA');
var kBldg    = makeKpi('BUILDINGS IN STUDY AREA');
var kExposed = makeKpi('IN ' + CFG.corridorM + ' m CORRIDOR');

body.add(ui.Panel([kReach.panel, kBasin.panel],
                  ui.Panel.Layout.flow('horizontal')));
body.add(ui.Panel([kBldg.panel, kExposed.panel],
                  ui.Panel.Layout.flow('horizontal')));

// ---- imagery controls ------------------------------------------------------
var imgCard = ui.Panel({style: {margin: '8px 4px 0 4px', padding: '9px',
                                backgroundColor: '#ffffff',
                                border: '1px solid ' + T.border}});
imgCard.add(ui.Label('Sentinel-2 imagery', {
  fontSize: '12px', fontWeight: 'bold', color: T.navy,
  backgroundColor: '#ffffff', margin: '0 0 6px 2px'
}));

var bandSelect = ui.Select({
  items: Object.keys(vizModes),
  value: defaultMode,
  style: {stretch: 'horizontal'}
});
imgCard.add(bandSelect);

var frameSelect = ui.Select({
  items: [{label: 'Cross-fade (use slider)',     value: 'blend'},
          {label: 'Post-event (Aug 27, 2026)',   value: 'post'},
          {label: 'Pre-event (1-month median)',  value: 'pre'},
          {label: 'Both (toggle in layer list)', value: 'both'},
          {label: 'Hide imagery',                value: 'none'}],
  value: 'blend',
  style: {stretch: 'horizontal', margin: '6px 0 0 0'}
});
imgCard.add(frameSelect);

// Cross-fade slider. layerPre sits UNDER layerPost in the stack, so lowering
// the post layer's opacity reveals the pre-event composite beneath it.
var blendRow = ui.Panel({style: {margin: '8px 0 0 0',
                                 backgroundColor: '#ffffff'}});
blendRow.add(ui.Label('Pre  ←—→  Post', {
  fontSize: '10px', color: T.textSub,
  backgroundColor: '#ffffff', margin: '0 0 2px 2px'
}));

var blendSlider = ui.Slider({
  min: 0, max: 1, value: 1, step: 0.02,
  style: {stretch: 'horizontal', margin: '0 2px'}
});
blendRow.add(blendSlider);
imgCard.add(blendRow);

var bandNote = ui.Label('', {fontSize: '10px', color: T.textSub,
                             backgroundColor: '#ffffff', margin: '6px 0 0 2px'});
imgCard.add(bandNote);
body.add(imgCard);

// ---- corridor control ------------------------------------------------------
var ctrl = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'),
                     style: {margin: '8px 0 0 4px', backgroundColor: T.cream}});
ctrl.add(ui.Label('Corridor width (m):', {fontSize: '12px', color: T.navy,
                                          backgroundColor: T.cream,
                                          margin: '8px 6px 0 0'}));
var corridorSelect = ui.Select({
  items: ['50', '100', '200', '300', '500'],
  value: String(CFG.corridorM), style: {width: '90px'}
});
ctrl.add(corridorSelect);
body.add(ctrl);

var status = ui.Label('Tracing downstream network…', {
  fontSize: '11px', color: T.teal, backgroundColor: T.cream, margin: '8px 6px'
});
body.add(status);

var chartSlot = ui.Panel({style: {margin: '4px 0 0 0',
                                  backgroundColor: T.cream}});
var stormSlot = ui.Panel({style: {margin: '0', backgroundColor: T.cream}});
body.add(chartSlot);
body.add(stormSlot);

var exportBtn = ui.Button({label: 'Queue exposure exports to Drive',
                           style: {stretch: 'horizontal', margin: '8px 4px'}});
body.add(exportBtn);

// ---- legend ----------------------------------------------------------------
var legend = ui.Panel({style: {margin: '8px 4px', padding: '8px',
                               backgroundColor: '#ffffff',
                               border: '1px solid ' + T.border}});
function legendRow(colour, text) {
  return ui.Panel([
    ui.Label('', {backgroundColor: colour, padding: '7px',
                  margin: '2px 6px 2px 2px', border: '1px solid ' + T.border}),
    ui.Label(text, {fontSize: '11px', color: T.navy,
                    backgroundColor: '#ffffff', margin: '3px 0'})
  ], ui.Panel.Layout.flow('horizontal'),
     {backgroundColor: '#ffffff', margin: '0'});
}

var bandLegend = ui.Panel({style: {backgroundColor: '#ffffff', margin: '0'}});

function renderBandLegend(key) {
  bandLegend.clear();
  bandLegend.add(ui.Label('Band → channel', {
    fontSize: '11px', fontWeight: 'bold', color: T.navy,
    backgroundColor: '#ffffff', margin: '2px 2px 5px 2px'
  }));
  var rows;
  if (key === 'swir') {
    rows = [['red', 'B11 (SWIR)'], ['green', 'B8 (NIR)'], ['blue', 'B4 (Red)']];
  } else if (key === 'fcc') {
    rows = [['red', 'B8 (NIR)'], ['green', 'B4 (Red)'], ['blue', 'B3 (Green)']];
  } else {
    rows = [['red', 'B4 (Red)'], ['green', 'B3 (Green)'], ['blue', 'B2 (Blue)']];
  }
  rows.forEach(function (r) { bandLegend.add(legendRow(r[0], r[1])); });
}

legend.add(bandLegend);
legend.add(ui.Panel({style: {height: '1px', backgroundColor: T.border,
                             margin: '7px 0', stretch: 'horizontal'}}));
legend.add(ui.Label('Map features', {fontSize: '11px', fontWeight: 'bold',
                                     color: T.navy, backgroundColor: '#ffffff',
                                     margin: '2px 2px 5px 2px'}));
legend.add(legendRow(T.teal,    'Downstream reaches (traced)'));
legend.add(legendRow(T.seafoam, 'Upstream / other reaches'));

// Corridor row lives in its own slot so the width can be re-rendered
var corridorLegendSlot = ui.Panel({style: {backgroundColor: '#ffffff',
                                           margin: '0'}});
legend.add(corridorLegendSlot);

function renderCorridorLegend() {
  corridorLegendSlot.clear();
  corridorLegendSlot.add(legendRow('#7fb3b3',
    CFG.corridorM + ' m downstream corridor'));
}

legend.add(legendRow(T.navy,  'Building footprint'));
legend.add(legendRow(T.alert, 'Footprint inside corridor'));
body.add(legend);

// ---- citation / provenance -------------------------------------------------
var citationBox = ui.Panel({style: {margin: '8px 4px 0 4px', padding: '9px',
                                    backgroundColor: '#ffffff',
                                    border: '1px solid ' + T.border}});

citationBox.add(ui.Label('Data sources', {
  fontSize: '11px', fontWeight: 'bold', color: T.navy,
  backgroundColor: '#ffffff', margin: '0 0 5px 2px'
}));
citationBox.add(ui.Label('Imagery: Copernicus Sentinel-2 MSI (S2_SR_HARMONIZED)', {
  fontSize: '10px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '0 0 2px 2px'
}));
citationBox.add(ui.Label('Cloud mask: Google Cloud Score+ (cs ≥ ' +
                         CFG.csThreshold + ')', {
  fontSize: '10px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '0 0 2px 2px'
}));
citationBox.add(ui.Label('Hydrology: HydroATLAS RiverATLAS v10 · HydroBASINS ' +
                         'level ' + CFG.basinLevel, {
  fontSize: '10px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '0 0 2px 2px'
}));
citationBox.add(ui.Label('Terrain: JAXA ALOS AW3D30 v4.1', {
  fontSize: '10px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '0 0 2px 2px'
}));
citationBox.add(ui.Label('Buildings: VIDA combined Google + Microsoft (' +
                         CFG.iso3 + ')', {
  fontSize: '10px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '0 0 2px 2px'
}));
citationBox.add(ui.Label('Rainfall: UCSB-CHG CHIRPS Pentad', {
  fontSize: '10px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '0 0 8px 2px'
}));

citationBox.add(ui.Panel({style: {height: '1px', backgroundColor: T.border,
                                  margin: '0 0 8px 0', stretch: 'horizontal'}}));

citationBox.add(ui.Label('Pennan Chinnasamy, Somdeep Kundu', {
  fontSize: '11px', fontWeight: 'bold', color: T.navy,
  backgroundColor: '#ffffff', margin: '0 0 2px 2px'
}));
citationBox.add(ui.Label('RuDRA Lab · C-TARA, IIT Bombay', {
  fontSize: '10px', fontWeight: 'bold', color: T.teal,
  backgroundColor: '#ffffff', margin: '0 0 0 2px'
}));

body.add(citationBox);

renderBandLegend('swir');
renderCorridorLegend();

// ============================================================================
// 7. MAP LAYERS
// ============================================================================
function drawLayers(corr) {
  Map.layers().reset([layerPre, layerPost]);   // imagery sits at the bottom

  Map.addLayer(basin.style({color: '457b9d', fillColor: '00000000', width: 2}),
    {}, 'Catchment', false);

  Map.addLayer(ee.FeatureCollection([ee.Feature(corr)])
      .style({color: '7fb3b3', fillColor: '2d8b8b40', width: 1}),
    {}, CFG.corridorM + ' m downstream corridor');

  Map.addLayer(upstream.style({color: 'a8dadc', width: 1}),
    {}, 'Upstream / other reaches', false);

  Map.addLayer(downstream.style({color: '2d8b8b', width: 2.5}),
    {}, 'Downstream reaches');

  Map.addLayer(buildings.style({color: '1a2332', fillColor: '00000000', width: 1}),
    {}, 'Buildings', false);

  Map.addLayer(exposed.style({color: 'e63946', fillColor: 'e6394688', width: 1}),
    {}, 'Buildings in corridor');

  Map.addLayer(ee.FeatureCollection([ee.Feature(aoi)])
      .style({color: 'f1faee', fillColor: '00000000', width: 2}), {}, 'Study area');

  Map.addLayer(ee.FeatureCollection([ee.Feature(pt)])
      .style({color: 'e63946', pointSize: 8, width: 2}), {}, 'Epicentre');
}

// ============================================================================
// 8. RAINFALL CHARTS
// ============================================================================
var MON = ['Jan','Feb','Mar','Apr','May','Jun',
           'Jul','Aug','Sep','Oct','Nov','Dec'];

function buildRainfallCharts(byYear) {
  chartSlot.clear();

  var head = ['Month'];
  CFG.years.forEach(function (y) { head.push(String(y)); });
  head.push('Mean ' + CFG.years[0] + '–' + (CFG.eventYear - 1));

  var rows = [head];
  for (var mi = 1; mi <= 12; mi++) {
    var key = 'c' + (mi < 10 ? '0' + mi : mi);
    var row = [MON[mi - 1]];
    var pSum = 0, pN = 0;

    CFG.years.forEach(function (y) {
      var d = byYear['y' + y];
      var v = (d && d[key] !== null && d[key] !== undefined) ? d[key] : null;
      if (y === CFG.eventYear && mi > CFG.monthsAvailable) { v = null; }
      row.push(v);
      if (y !== CFG.eventYear && v !== null) { pSum += v; pN++; }
    });
    row.push(pN ? pSum / pN : null);
    rows.push(row);
  }

  var colors = CFG.years.map(function (y) {
    return y === CFG.eventYear ? T.alert : T.seafoam;
  });
  colors.push(T.navy);

  var series = {};
  CFG.years.forEach(function (y, i) {
    series[i] = (y === CFG.eventYear)
      ? {lineWidth: 3, pointSize: 5} : {lineWidth: 1, pointSize: 0};
  });
  series[CFG.years.length] = {lineWidth: 2, lineDashStyle: [6, 4], pointSize: 0};

  chartSlot.add(ui.Chart(rows, 'LineChart', {
    title: 'Cumulative rainfall by year — catchment mean (CHIRPS)',
    titleTextStyle: {fontSize: 12, color: T.navy, bold: true},
    hAxis: {title: 'Month', textStyle: {color: T.textSub}},
    vAxis: {title: 'Cumulative rainfall (mm)', textStyle: {color: T.textSub}},
    colors: colors, series: series, backgroundColor: T.cream,
    legend: {position: 'right', textStyle: {fontSize: 10, color: T.navy}},
    chartArea: {width: '66%', height: '70%'},
    height: 260
  }));

  var combo = [['Month', CFG.eventYear + ' monthly',
                CFG.eventYear + ' cumulative', 'Normal cumulative']];
  for (var m2 = 1; m2 <= 12; m2++) {
    var mk = 'm' + (m2 < 10 ? '0' + m2 : m2);
    var ck = 'c' + (m2 < 10 ? '0' + m2 : m2);
    var ev = byYear['y' + CFG.eventYear] || {};
    var monthly = (m2 <= CFG.monthsAvailable) ? (ev[mk] || null) : null;
    var cumul   = (m2 <= CFG.monthsAvailable) ? (ev[ck] || null) : null;

    var nSum = 0, nN = 0;
    CFG.years.forEach(function (y) {
      if (y === CFG.eventYear) { return; }
      var d = byYear['y' + y];
      if (d && d[ck] !== null && d[ck] !== undefined) { nSum += d[ck]; nN++; }
    });
    combo.push([MON[m2 - 1], monthly, cumul, nN ? nSum / nN : null]);
  }

  chartSlot.add(ui.Chart(combo, 'ComboChart', {
    title: CFG.eventYear + ' against the ' + (CFG.years.length - 1) + '-year normal',
    titleTextStyle: {fontSize: 12, color: T.navy, bold: true},
    hAxis: {title: 'Month', textStyle: {color: T.textSub}},
    vAxes: {0: {title: 'Monthly (mm)', gridlines: {color: 'none'}},
            1: {title: 'Cumulative (mm)'}},
    series: {
      0: {type: 'bars', color: T.teal, targetAxisIndex: 0},
      1: {type: 'line', color: T.alert, lineWidth: 3, pointSize: 4,
          targetAxisIndex: 1},
      2: {type: 'line', color: T.navy, lineWidth: 2, pointSize: 0,
          lineDashStyle: [6, 4], targetAxisIndex: 1}
    },
    backgroundColor: T.cream,
    legend: {position: 'top', alignment: 'end',
             textStyle: {fontSize: 10, color: T.navy}},
    chartArea: {width: '66%', height: '60%'},
    height: 250
  }));
}

function buildStormChart() {
  stormSlot.clear();
  stormSlot.add(ui.Chart.image.series({
    imageCollection: stormColl, region: basinGeom,
    reducer: ee.Reducer.mean(), scale: CFG.chirpsScale,
    xProperty: 'system:time_start'
  }).setChartType('ColumnChart').setOptions({
    title: 'Pentad rainfall around the event',
    titleTextStyle: {fontSize: 12, color: T.navy, bold: true},
    hAxis: {title: 'Date', format: 'dd MMM', textStyle: {color: T.textSub}},
    vAxis: {title: 'Rainfall (mm / 5 days)', textStyle: {color: T.textSub}},
    colors: [T.teal], backgroundColor: T.cream,
    legend: {position: 'none'},
    chartArea: {width: '74%', height: '60%'},
    height: 200
  }));
}

// ============================================================================
// 9. MAIN
// ============================================================================
var pending = null;

function run() {
  status.setValue('Computing…');
  status.style().set('color', T.teal);
  chartSlot.clear();
  stormSlot.clear();

  // reflect the current corridor width in the KPI title and the legend
  kExposed.title.setValue('IN ' + CFG.corridorM + ' m CORRIDOR');
  renderCorridorLegend();

  corridor = downstream.geometry(10)
               .buffer(CFG.corridorM, 10)
               .intersection(aoi, 10);
  exposed = buildings.filterBounds(corridor);

  drawLayers(corridor);
  applyFrame(frameSelect.getValue());

  ee.Dictionary({
    reachKm:  downstream.geometry(10).length(10).divide(1000),
    nReach:   downstream.size(),
    basinKm2: basinGeom.area(30).divide(1e6),
    nBldg:    buildings.size(),
    nExposed: exposed.size(),
    roofExp:  exposed.aggregate_sum('area_in_meters'),
    drop:     seedElev.subtract(ee.Number(downstream.aggregate_min('_elev'))),
    nPre:     preCollection.size()
  }).evaluate(function (d, err) {
    if (err || !d) {
      status.setValue('Stats failed: ' + (err || 'no result'));
      status.style().set('color', T.alert);
      return;
    }
    var f = function (x, dp) {
      return (x === null || x === undefined) ? 'n/a' :
        Number(x).toFixed(dp === undefined ? 0 : dp)
          .replace(/\B(?=(\d{3})+(?!\d))/g, ',');
    };

    kReach.value.setValue(f(d.reachKm, 1) + ' km');
    kReach.sub.setValue(f(d.nReach) + ' reaches · ' + f(d.drop) + ' m descent');

    kBasin.value.setValue(f(d.basinKm2, 1) + ' km²');
    kBasin.sub.setValue('HydroBASINS level ' + CFG.basinLevel);

    kBldg.value.setValue(f(d.nBldg));
    kBldg.sub.setValue('VIDA ' + CFG.iso3 + ' · ' +
                       (CFG.aoiRadiusM / 1000) + ' km study area');

    kExposed.value.setValue(f(d.nExposed));
    var share = d.nBldg ? (100 * d.nExposed / d.nBldg) : 0;
    kExposed.sub.setValue(share.toFixed(1) + '% of stock · ' +
                          f((d.roofExp || 0) / 1000, 1) + 'k m² roof');
    kExposed.value.style().set('color', share > 15 ? T.alert : T.teal);

    bandNote.setValue('Pre-event composite built from ' + f(d.nPre) +
                      ' scenes · post-event 2026-08-27');
  });

  var bundle = ee.Dictionary({});
  CFG.years.forEach(function (y) { bundle = bundle.set('y' + y, yearStats(y)); });

  bundle.evaluate(function (d, err) {
    if (err || !d) {
      status.setValue('Rainfall failed: ' + (err || 'no result'));
      status.style().set('color', T.alert);
      return;
    }
    buildRainfallCharts(d);
    buildStormChart();
    status.setValue('Done · downstream trace within the study area, CHIRPS ' +
                    CFG.years[0] + '–' + CFG.eventYear);
    status.style().set('color', T.teal);
  });

  pending = {downstream: downstream, corridor: corridor, exposed: exposed,
             basin: basin};
}

// ============================================================================
// 10. INTERACTION
// ============================================================================
function applyFrame(mode) {
  if (mode === 'blend') {
    layerPre.setShown(true);
    layerPost.setShown(true);
    layerPre.setOpacity(1);
    layerPost.setOpacity(blendSlider.getValue());
    blendRow.style().set('shown', true);
  } else {
    layerPre.setOpacity(1);
    layerPost.setOpacity(1);
    layerPre.setShown(mode === 'pre' || mode === 'both');
    layerPost.setShown(mode === 'post' || mode === 'both');
    blendRow.style().set('shown', false);
  }
}

blendSlider.onChange(function (v) {
  if (frameSelect.getValue() === 'blend') { layerPost.setOpacity(v); }
});

bandSelect.onChange(function (k) {
  var v = vizModes[k];
  layerPre.setVisParams(v);
  layerPost.setVisParams(v);
  renderBandLegend(v.key);
});

frameSelect.onChange(applyFrame);

corridorSelect.onChange(function (v) {
  CFG.corridorM = parseInt(v, 10);
  kExposed.value.setValue('…');
  kExposed.sub.setValue(' ');
  run();
});

exportBtn.onClick(function () {
  if (!pending) { return; }
  var tag = 'rasuwa_ds_' + CFG.corridorM + 'm';

  Export.table.toDrive({collection: pending.exposed,
    description: 'exposed_buildings_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});
  Export.table.toDrive({collection: pending.downstream,
    description: 'downstream_reaches_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});
  Export.table.toDrive({
    collection: ee.FeatureCollection([ee.Feature(pending.corridor)]),
    description: 'downstream_corridor_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});

  status.setValue('3 tasks queued — open the Tasks tab and press Run on each.');
  status.style().set('color', T.navy);
});

// ============================================================================
// 11. GO
// ============================================================================
print('Seed reach (should sit on the epicentre):', seed);
print('Downstream reach IDs (full trace):', chain);
print('Downstream reaches kept in study area:', downstream.size());
print('Pre-event scenes in composite window:', preCollection.size());
print('Building sample (verify VIDA ' + CFG.iso3 + ' resolved):', buildings.first());

run();
