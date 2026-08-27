/*******************************************************************************
 * NEPAL DEBRIS FLOW — DOWNSTREAM EXPOSURE, IMAGERY & RAINFALL ANOMALY  (v3.2)
 * Rasuwa River Valley (28.279 N, 85.378 E) — event 2026-08-25
 *
 * v3.2 optimizations — reduced hanging/timeouts on heavy computation:
 *   • vectorization: maxPixels 1e10 → 1e8, tileScale 4 → 2
 *   • HAND inundation: tileScale 4 → 2, scale 90 → 100, eightConnected off
 *   • SAR filtering: clipped to AOI before averaging
 *
 * v3 adds, on top of the original downstream-exposure browser:
 *   1. PlanetScope high-resolution pre/post mosaics (own frame + band control)
 *   2. Change detection — Sentinel-2 dNDVI disturbance footprint (cloud-masked)
 *                       — Sentinel-1 VV backscatter change (cloud-proof SAR)
 *   3. Broadened exposure — WorldPop population + ESA WorldCover cropland
 *   4. Terrain-aware hazard — MERIT Hydro HAND inundation, selectable as the
 *      exposure footprint (River buffer / HAND / Observed disturbance)
 *   5. Storm intensity — GPM IMERG pre-event rainfall total
 *   6. Field-validation & glacial-lake hooks (paste an asset ID to enable)
 *   7. Expanded KPIs, legend, exports and provenance
 *
 * Imagery   : S2 pre 1-month median vs post scene; PlanetScope ~3 m pre/post
 * Hydrology : downstream-only trace via RiverATLAS NEXT_DOWN, clipped to AOI
 * Rainfall  : CHIRPS cumulative (context) + GPM IMERG (event intensity)
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
  textSub: '#5c6b73',
  // v3 accents
  amber:   '#e9c46a',   // observed disturbance
  handblue:'#457b9d',   // HAND inundation
  sar:     '#8d5bd6',   // SAR change
  pop:     '#e76f51'    // population
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
  monthsAvailable: 8,

  // ---- v3 additions --------------------------------------------------------
  hazardBasis: 'buffer',   // 'buffer' | 'hand' | 'observed'  (default preserves v2)

  // Change detection (Sentinel-2, both surface reflectance → comparable)
  dndviDrop:   -0.15,      // NDVI drop that flags vegetation loss / scour
  briRise:      600,       // SWIR/red brightness rise that flags fresh debris (B11)
  changeScale:  20,

  // Sentinel-1 SAR
  s1PreDays:  24,          // pre-event window for the SAR baseline
  s1PostDays: 8,           // post-event window
  s1ChangeDb: -3,          // VV drop (dB) that flags smooth/wet new surface

  // Terrain-aware inundation (MERIT Hydro HAND)
  handThreshM: 15,         // height-above-drainage cut for potential inundation
  handBufM:    1500,       // confine HAND zone to this buffer of the traced channel

  // Population & land cover
  popYear:  2020,          // WorldPop GP 100 m latest reliable year

  // Storm intensity (GPM IMERG)
  imergPreDays: 3,         // hours-before window for the triggering storm total

  // Optional user assets — paste an EE asset ID to switch these on
  fieldPointsAsset: '',    // e.g. 'users/you/rasuwa_field_points' (geotagged photos)
  glacialLakeAsset: ''     // e.g. 'users/you/icimod_glacial_lakes'
};

var EVENT = ee.Date(CFG.eventDate);
var pt = ee.Geometry.Point([CFG.lon, CFG.lat]);
var aoi = pt.buffer(CFG.aoiRadiusM);
var searchArea = pt.buffer(CFG.searchM);

// ============================================================================
// 2. IMAGERY — Sentinel-2 pre-event composite & post-event scene
// ============================================================================
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var preCollection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(EVENT.advance(-CFG.baselineDays, 'day'), EVENT)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80));

var preMasked = preCollection.linkCollection(csPlus, ['cs'])
  .map(function (img) {
    return img.updateMask(img.select('cs').gte(CFG.csThreshold));
  });

var preImage  = preMasked.median().clip(aoi);
var postImage = ee.Image(CFG.postAssetId).clip(aoi);

// Cloud-masked post scene, used ONLY for change detection (display keeps the
// full scene). Link the single post image to Cloud Score+ by system:index.
var postMaskedImg = ee.Image(
  ee.ImageCollection([ee.Image(CFG.postAssetId)])
    .linkCollection(csPlus, ['cs'])
    .first());
var postCloudMask = postMaskedImg.select('cs').gte(CFG.csThreshold);
var postForChange = postMaskedImg.updateMask(postCloudMask).clip(aoi);

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
                             'Pre-event composite (S2)', false);
var layerPost = ui.Map.Layer(postImage, vizModes[defaultMode],
                             'Post-event S2 (Aug 27)', true);

// ============================================================================
// 2b. PLANETSCOPE — very-high-resolution pre/post mosaics (~3 m)
// ----------------------------------------------------------------------------
// NOTE ON RADIOMETRY: the pre tiles are *analytic_sr* (surface reflectance,
// 0–10000) while the post tiles are *analytic* (at-sensor radiance, a different
// scale). They are therefore NOT radiometrically comparable — Planet is used
// here as a high-resolution VISUAL layer only. Quantitative change detection
// uses Sentinel-2 (both SR) and Sentinel-1 (SAR). Each Planet epoch gets its
// own stretch below; tune the min/max to taste with the layer's own settings.
// PlanetScope 4-band order: b1 Blue · b2 Green · b3 Red · b4 NIR.
// ============================================================================
var planetPre = ee.ImageCollection([
  ee.Image('projects/ee-arpandawn/assets/20260527_053226_41_254a_analytic_sr'),
  ee.Image('projects/ee-arpandawn/assets/20260527_053224_18_254a_analytic_sr'),
  ee.Image('projects/ee-arpandawn/assets/20260527_053221_96_254a_analytic_sr'),
  ee.Image('projects/ee-arpandawn/assets/20260527_053219_95_254a_analytic_sr'),
  ee.Image('projects/ee-arpandawn/assets/20260527_053217_72_254a_analytic_sr')
]).mosaic();

var planetPost = ee.ImageCollection([
  ee.Image('projects/ee-arpandawn/assets/20260826_054502_86_251f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_054500_80_251f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_054458_74_251f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_054456_67_251f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_050135_34_255f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_050133_00_255f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_050130_66_255f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_050128_33_255f_analytic'),
  ee.Image('projects/ee-arpandawn/assets/20260826_050125_99_255f_analytic')
]).mosaic();

// Separate stretches: pre = SR, post = radiance.
var planetViz = {
  natural: {
    pre:  {bands: ['b3', 'b2', 'b1'], min: 200, max: 2500, gamma: 1.3},
    post: {bands: ['b3', 'b2', 'b1'], min: 300, max: 6000, gamma: 1.3}
  },
  cir: {   // colour-infrared: NIR · Red · Green — vegetation reads red
    pre:  {bands: ['b4', 'b3', 'b2'], min: 200, max: 4000, gamma: 1.2},
    post: {bands: ['b4', 'b3', 'b2'], min: 328, max: 9591, gamma: 1.0}
  }
};
var planetBand = 'natural';

var layerPlanetPre  = ui.Map.Layer(planetPre,  planetViz.natural.pre,
                                   'PlanetScope pre (May 27)', false, 0);
var layerPlanetPost = ui.Map.Layer(planetPost, planetViz.natural.post,
                                   'PlanetScope post (Aug 26)', false, 0);

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
var slope = ee.Terrain.slope(dem);   // v3: reused for the susceptibility hint

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
// 4b. CHANGE DETECTION  (optical dNDVI + SAR VV)  → observed hazard footprint
// ============================================================================
// --- Sentinel-2 vegetation-loss / scour signal (both scenes are SR) ---------
var preNDVI  = preImage.normalizedDifference(['B8', 'B4']).rename('ndvi');
var postNDVI = postForChange.normalizedDifference(['B8', 'B4']).rename('ndvi');
var dNDVI    = postNDVI.subtract(preNDVI).rename('dndvi');

// fresh debris also brightens in SWIR; combine the two evidence rules
var dBRI = postForChange.select('B11').subtract(preImage.select('B11')).rename('dbri');
var disturbedMask = dNDVI.lt(CFG.dndviDrop)
                      .and(slope.lt(35))            // drop steep-cloud artefacts
                      .rename('dist');
var disturbed = disturbedMask.selfMask();
var disturbAreaImg = disturbedMask.multiply(ee.Image.pixelArea()).rename('darea');

// vectorise the disturbance so it can drive exposure and be exported
// OPTIMIZED: reduce maxPixels, tileScale for faster processing
var disturbVec = disturbedMask.selfMask().reduceToVectors({
  geometry: aoi, scale: CFG.changeScale, geometryType: 'polygon',
  eightConnected: true, maxPixels: 1e8, tileScale: 2
});

// --- Sentinel-1 SAR backscatter change (cloud-proof) ------------------------
// OPTIMIZED: clip early to reduce processing footprint
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(aoi)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV');

var s1Pre  = s1.filterDate(EVENT.advance(-CFG.s1PreDays, 'day'), EVENT).median().clip(aoi);
var s1Post = s1.filterDate(EVENT, EVENT.advance(CFG.s1PostDays, 'day')).median().clip(aoi);
var s1Diff = s1Post.subtract(s1Pre).rename('vv_db');   // dB change
var s1NewSmooth = s1Diff.lt(CFG.s1ChangeDb).selfMask();          // wet / smooth

// ============================================================================
// 4c. TERRAIN-AWARE INUNDATION  (MERIT Hydro HAND)
// ============================================================================
var hand = ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var handMask = hand.lte(CFG.handThreshM);
// OPTIMIZED: reduced scale (100m→90m is fine), tileScale 2, maxPixels 5e7
var handZone = handMask.selfMask().reduceToVectors({
  geometry: downstream.geometry(10).buffer(CFG.handBufM, 10).intersection(aoi, 10),
  scale: 100, geometryType: 'polygon', eightConnected: false,
  maxPixels: 5e7, tileScale: 2
});
var handZoneGeom = handZone.geometry();

// ============================================================================
// 4d. BROADENED EXPOSURE  (population + cropland)
// ============================================================================
var pop = ee.ImageCollection('WorldPop/GP/100m/pop')
  .filter(ee.Filter.eq('year', CFG.popYear))
  .filterBounds(aoi)
  .mosaic()
  .rename('pop');

var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');
var cropland   = worldcover.eq(40).rename('crop');            // class 40 = cropland
var cropAreaImg = cropland.multiply(ee.Image.pixelArea()).rename('carea');

// ============================================================================
// 4e. OPTIONAL USER LAYERS  (field photos, glacial lakes)
// ============================================================================
var fieldPoints = CFG.fieldPointsAsset ?
  ee.FeatureCollection(CFG.fieldPointsAsset) : null;
var glacialLakes = CFG.glacialLakeAsset ?
  ee.FeatureCollection(CFG.glacialLakeAsset) : null;

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

// ---- GPM IMERG — event-window intensity (finer than CHIRPS pentad) ---------
// V07 half-hourly precipitation is in mm/hr; multiply each step by 0.5 h.
var imerg = ee.ImageCollection('NASA/GPM_L3/IMERG_V07')
  .filterDate(EVENT.advance(-CFG.imergPreDays, 'day'), EVENT.advance(1, 'day'))
  .select('precipitation');
var imergTotal = imerg.sum().multiply(0.5).rename('mm').clip(aoi);
var imergScale = 11132;

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
var kPeople  = makeKpi('PEOPLE IN HAZARD ZONE');     // v3
var kCrop    = makeKpi('CROPLAND IN HAZARD ZONE');   // v3
var kDisturb = makeKpi('OBSERVED DISTURBANCE');       // v3
var kRain    = makeKpi(CFG.imergPreDays + '-DAY STORM TOTAL'); // v3

body.add(ui.Panel([kReach.panel, kBasin.panel],
                  ui.Panel.Layout.flow('horizontal')));
body.add(ui.Panel([kBldg.panel, kExposed.panel],
                  ui.Panel.Layout.flow('horizontal')));
body.add(ui.Panel([kPeople.panel, kCrop.panel],
                  ui.Panel.Layout.flow('horizontal')));
body.add(ui.Panel([kDisturb.panel, kRain.panel],
                  ui.Panel.Layout.flow('horizontal')));

// ---- imagery controls (Sentinel-2) -----------------------------------------
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

// ---- PlanetScope high-res controls (v3) ------------------------------------
var planetCard = ui.Panel({style: {margin: '8px 4px 0 4px', padding: '9px',
                                   backgroundColor: '#ffffff',
                                   border: '1px solid ' + T.border}});
planetCard.add(ui.Label('PlanetScope ~3 m imagery', {
  fontSize: '12px', fontWeight: 'bold', color: T.navy,
  backgroundColor: '#ffffff', margin: '0 0 6px 2px'
}));

var planetBandSelect = ui.Select({
  items: [{label: 'Natural colour [b3,b2,b1]', value: 'natural'},
          {label: 'Colour-infrared [b4,b3,b2]', value: 'cir'}],
  value: 'natural', style: {stretch: 'horizontal'}
});
planetCard.add(planetBandSelect);

var planetFrameSelect = ui.Select({
  items: [{label: 'Hide PlanetScope',        value: 'none'},
          {label: 'Post (Aug 26, 2026)',     value: 'post'},
          {label: 'Pre (May 27, 2026)',      value: 'pre'},
          {label: 'Both (toggle in layers)', value: 'both'}],
  value: 'none', style: {stretch: 'horizontal', margin: '6px 0 0 0'}
});
planetCard.add(planetFrameSelect);
planetCard.add(ui.Label('Pre = surface reflectance, post = radiance — a visual '
  + 'reference, not radiometrically matched. Change stats use S2 & S1.', {
  fontSize: '9.5px', color: T.textSub, backgroundColor: '#ffffff',
  margin: '6px 0 0 2px'
}));
body.add(planetCard);

// ---- change-detection controls (v3) ----------------------------------------
var changeCard = ui.Panel({style: {margin: '8px 4px 0 4px', padding: '9px',
                                   backgroundColor: '#ffffff',
                                   border: '1px solid ' + T.border}});
changeCard.add(ui.Label('Change detection', {
  fontSize: '12px', fontWeight: 'bold', color: T.navy,
  backgroundColor: '#ffffff', margin: '0 0 6px 2px'
}));
var chkDndvi = ui.Checkbox({label: 'S2 disturbance (dNDVI, cloud-sensitive)',
  value: true, style: {fontSize: '11px', backgroundColor: '#ffffff'}});
var chkSar   = ui.Checkbox({label: 'S1 SAR change (VV, cloud-proof)',
  value: true, style: {fontSize: '11px', backgroundColor: '#ffffff'}});
var chkHand  = ui.Checkbox({label: 'HAND inundation (≤ ' + CFG.handThreshM + ' m)',
  value: false, style: {fontSize: '11px', backgroundColor: '#ffffff'}});
changeCard.add(chkDndvi);
changeCard.add(chkSar);
changeCard.add(chkHand);
body.add(changeCard);

// ---- hazard-basis + corridor control ---------------------------------------
var ctrl = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'),
                     style: {margin: '8px 0 0 4px', backgroundColor: T.cream}});
ctrl.add(ui.Label('Exposure footprint:', {fontSize: '12px', color: T.navy,
                                          backgroundColor: T.cream,
                                          margin: '8px 6px 0 0'}));
var hazardSelect = ui.Select({
  items: [{label: 'River buffer',        value: 'buffer'},
          {label: 'HAND inundation',     value: 'hand'},
          {label: 'Observed disturbance',value: 'observed'}],
  value: CFG.hazardBasis, style: {width: '150px'}
});
ctrl.add(hazardSelect);
body.add(ctrl);

var ctrl2 = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'),
                      style: {margin: '6px 0 0 4px', backgroundColor: T.cream}});
ctrl2.add(ui.Label('Corridor width (m):', {fontSize: '12px', color: T.navy,
                                          backgroundColor: T.cream,
                                          margin: '8px 6px 0 0'}));
var corridorSelect = ui.Select({
  items: ['50', '100', '200', '300', '500'],
  value: String(CFG.corridorM), style: {width: '90px'}
});
ctrl2.add(corridorSelect);
body.add(ctrl2);

var status = ui.Label('Tracing downstream network…', {
  fontSize: '11px', color: T.teal, backgroundColor: T.cream, margin: '8px 6px'
});
body.add(status);

var chartSlot = ui.Panel({style: {margin: '4px 0 0 0',
                                  backgroundColor: T.cream}});
var stormSlot = ui.Panel({style: {margin: '0', backgroundColor: T.cream}});
body.add(chartSlot);
body.add(stormSlot);

var exportBtn = ui.Button({label: 'Queue exposure + hazard exports to Drive',
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
  bandLegend.add(ui.Label('S2 band → channel', {
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

legend.add(legendRow(T.navy,     'Building footprint'));
legend.add(legendRow(T.alert,    'Footprint inside hazard zone'));
legend.add(legendRow(T.amber,    'Observed disturbance (S2 dNDVI)'));
legend.add(legendRow(T.sar,      'SAR change — wet/smooth (S1 VV)'));
legend.add(legendRow(T.handblue, 'HAND inundation (≤ ' + CFG.handThreshM + ' m)'));
body.add(legend);

// ---- citation / provenance -------------------------------------------------
var citationBox = ui.Panel({style: {margin: '8px 4px 0 4px', padding: '9px',
                                    backgroundColor: '#ffffff',
                                    border: '1px solid ' + T.border}});

citationBox.add(ui.Label('Data sources', {
  fontSize: '11px', fontWeight: 'bold', color: T.navy,
  backgroundColor: '#ffffff', margin: '0 0 5px 2px'
}));
function srcLine(txt) {
  citationBox.add(ui.Label(txt, {fontSize: '10px', color: T.textSub,
    backgroundColor: '#ffffff', margin: '0 0 2px 2px'}));
}
srcLine('Imagery: Copernicus Sentinel-2 MSI (S2_SR_HARMONIZED)');
srcLine('High-res: PlanetScope analytic / SR (© Planet Labs PBC)');
srcLine('SAR: Copernicus Sentinel-1 GRD (VV)');
srcLine('Cloud mask: Google Cloud Score+ (cs ≥ ' + CFG.csThreshold + ')');
srcLine('Hydrology: HydroATLAS RiverATLAS v10 · HydroBASINS lvl ' + CFG.basinLevel);
srcLine('Inundation: MERIT Hydro v1.0.1 (HAND)');
srcLine('Terrain: JAXA ALOS AW3D30 v4.1');
srcLine('Buildings: VIDA combined Google + Microsoft (' + CFG.iso3 + ')');
srcLine('Population: WorldPop GP 100 m (' + CFG.popYear + ')');
srcLine('Land cover: ESA WorldCover v200');
srcLine('Rainfall: UCSB-CHG CHIRPS Pentad · NASA GPM IMERG V07');

citationBox.add(ui.Panel({style: {height: '1px', backgroundColor: T.border,
                                  margin: '6px 0 8px 0', stretch: 'horizontal'}}));

citationBox.add(ui.Label('Pennan Chinnasamy, Surajit Ghosh, Somdeep Kundu', {
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
function drawLayers(hazardGeom) {
  // imagery sits at the bottom: Planet under S2
  Map.layers().reset([layerPlanetPre, layerPlanetPost, layerPre, layerPost]);

  Map.addLayer(basin.style({color: '457b9d', fillColor: '00000000', width: 2}),
    {}, 'Catchment', false);

  if (chkHand.getValue()) {
    Map.addLayer(handMask.selfMask(),
      {palette: [T.handblue.replace('#', '')]}, 'HAND inundation', true, 0.5);
  }

  Map.addLayer(ee.FeatureCollection([ee.Feature(hazardGeom)])
      .style({color: '7fb3b3', fillColor: '2d8b8b40', width: 1}),
    {}, 'Hazard footprint (' + CFG.hazardBasis + ')');

  Map.addLayer(upstream.style({color: 'a8dadc', width: 1}),
    {}, 'Upstream / other reaches', false);

  Map.addLayer(downstream.style({color: '2d8b8b', width: 2.5}),
    {}, 'Downstream reaches');

  if (chkSar.getValue()) {
    Map.addLayer(s1Diff, {min: -6, max: 6,
      palette: ['8d5bd6', 'ffffff', '2d8b8b']}, 'SAR VV change (dB)', false);
    Map.addLayer(s1NewSmooth, {palette: [T.sar.replace('#', '')]},
      'SAR: new wet/smooth', true, 0.7);
  }

  if (chkDndvi.getValue()) {
    Map.addLayer(disturbed, {palette: [T.amber.replace('#', '')]},
      'Observed disturbance (dNDVI)', true, 0.75);
  }

  Map.addLayer(buildings.style({color: '1a2332', fillColor: '00000000', width: 1}),
    {}, 'Buildings', false);

  Map.addLayer(exposed.style({color: 'e63946', fillColor: 'e6394688', width: 1}),
    {}, 'Buildings in hazard zone');

  if (fieldPoints) {
    Map.addLayer(fieldPoints.style({color: 'ffd166', pointSize: 6, width: 2}),
      {}, 'Field validation points');
  }
  if (glacialLakes) {
    Map.addLayer(glacialLakes.style({color: '48cae4', fillColor: '48cae466',
      width: 1}), {}, 'Glacial lakes', false);
  }

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
    title: 'Pentad rainfall around the event (CHIRPS)',
    titleTextStyle: {fontSize: 12, color: T.navy, bold: true},
    hAxis: {title: 'Date', format: 'dd MMM', textStyle: {color: T.textSub}},
    vAxis: {title: 'Rainfall (mm / 5 days)', textStyle: {color: T.textSub}},
    colors: [T.teal], backgroundColor: T.cream,
    legend: {position: 'none'},
    chartArea: {width: '74%', height: '60%'},
    height: 200
  }));

  // v3: IMERG half-hourly intensity around the event
  stormSlot.add(ui.Chart.image.series({
    imageCollection: ee.ImageCollection('NASA/GPM_L3/IMERG_V07')
      .filterDate(EVENT.advance(-7, 'day'), EVENT.advance(3, 'day'))
      .select('precipitation'),
    region: basinGeom, reducer: ee.Reducer.mean(), scale: imergScale,
    xProperty: 'system:time_start'
  }).setChartType('ColumnChart').setOptions({
    title: 'GPM IMERG half-hourly intensity around the event',
    titleTextStyle: {fontSize: 12, color: T.navy, bold: true},
    hAxis: {title: 'Date', format: 'dd MMM', textStyle: {color: T.textSub}},
    vAxis: {title: 'Rain rate (mm / hr)', textStyle: {color: T.textSub}},
    colors: [T.alert], backgroundColor: T.cream,
    legend: {position: 'none'},
    chartArea: {width: '74%', height: '60%'},
    height: 200
  }));
}

// ============================================================================
// 9. MAIN
// ============================================================================
var pending = null;

// Resolve the exposure footprint geometry from the selected basis
function resolveHazardGeom() {
  if (CFG.hazardBasis === 'observed') { return disturbVec.geometry(); }
  if (CFG.hazardBasis === 'hand')     { return handZoneGeom; }
  return downstream.geometry(10).buffer(CFG.corridorM, 10).intersection(aoi, 10);
}

function run() {
  status.setValue('Computing…');
  status.style().set('color', T.teal);
  chartSlot.clear();
  stormSlot.clear();

  CFG.hazardBasis = hazardSelect.getValue();

  // reflect current settings in KPI titles and legend
  kExposed.title.setValue(CFG.hazardBasis === 'buffer'
    ? 'IN ' + CFG.corridorM + ' m CORRIDOR'
    : 'IN HAZARD ZONE (' + CFG.hazardBasis + ')');
  renderCorridorLegend();

  var hazardGeom = resolveHazardGeom();
  corridor = hazardGeom;
  exposed = buildings.filterBounds(hazardGeom);

  drawLayers(hazardGeom);
  applyFrame(frameSelect.getValue());
  applyPlanetFrame(planetFrameSelect.getValue());

  ee.Dictionary({
    reachKm:   downstream.geometry(10).length(10).divide(1000),
    nReach:    downstream.size(),
    basinKm2:  basinGeom.area(30).divide(1e6),
    nBldg:     buildings.size(),
    nExposed:  exposed.size(),
    roofExp:   exposed.aggregate_sum('area_in_meters'),
    drop:      seedElev.subtract(ee.Number(downstream.aggregate_min('_elev'))),
    nPre:      preCollection.size(),
    // v3 metrics
    people:    pop.reduceRegion({reducer: ee.Reducer.sum(), geometry: hazardGeom,
                 scale: 100, maxPixels: 1e9, bestEffort: true}).get('pop'),
    cropHa:    ee.Number(cropAreaImg.reduceRegion({reducer: ee.Reducer.sum(),
                 geometry: hazardGeom, scale: 10, maxPixels: 1e10,
                 bestEffort: true, tileScale: 4}).get('carea')).divide(1e4),
    distHa:    ee.Number(disturbAreaImg.reduceRegion({reducer: ee.Reducer.sum(),
                 geometry: aoi, scale: CFG.changeScale, maxPixels: 1e10,
                 bestEffort: true, tileScale: 4}).get('darea')).divide(1e4),
    rainMm:    imergTotal.reduceRegion({reducer: ee.Reducer.mean(),
                 geometry: basinGeom, scale: imergScale, maxPixels: 1e9,
                 bestEffort: true}).get('mm')
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

    kPeople.value.setValue(f(d.people, 0));
    kPeople.sub.setValue('WorldPop ' + CFG.popYear + ' · ' + CFG.hazardBasis);
    kPeople.value.style().set('color', (d.people || 0) > 500 ? T.alert : T.teal);

    kCrop.value.setValue(f(d.cropHa, 1) + ' ha');
    kCrop.sub.setValue('ESA WorldCover cropland');

    kDisturb.value.setValue(f(d.distHa, 1) + ' ha');
    kDisturb.sub.setValue('S2 dNDVI < ' + CFG.dndviDrop + ' (cloud-sensitive)');

    kRain.value.setValue(f(d.rainMm, 0) + ' mm');
    kRain.sub.setValue('IMERG ' + CFG.imergPreDays + '-day pre-event · basin mean');
    kRain.value.style().set('color', (d.rainMm || 0) > 150 ? T.alert : T.teal);

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
    status.setValue('Done · ' + CFG.hazardBasis + ' footprint · CHIRPS ' +
                    CFG.years[0] + '–' + CFG.eventYear + ' + IMERG');
    status.style().set('color', T.teal);
  });

  pending = {downstream: downstream, hazardGeom: hazardGeom, exposed: exposed,
             basin: basin, disturbVec: disturbVec, handZone: handZone};
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

function applyPlanetFrame(mode) {
  layerPlanetPre.setShown(mode === 'pre' || mode === 'both');
  layerPlanetPost.setShown(mode === 'post' || mode === 'both');
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

planetBandSelect.onChange(function (k) {
  planetBand = k;
  layerPlanetPre.setVisParams(planetViz[k].pre);
  layerPlanetPost.setVisParams(planetViz[k].post);
});
planetFrameSelect.onChange(applyPlanetFrame);

hazardSelect.onChange(function () {
  kPeople.value.setValue('…'); kCrop.value.setValue('…');
  run();
});

corridorSelect.onChange(function (v) {
  CFG.corridorM = parseInt(v, 10);
  kExposed.value.setValue('…');
  kExposed.sub.setValue(' ');
  run();
});

// change-detection toggles just re-draw the map (stats are already computed)
chkDndvi.onChange(function () { drawLayers(resolveHazardGeom());
  applyFrame(frameSelect.getValue()); applyPlanetFrame(planetFrameSelect.getValue()); });
chkSar.onChange(function ()   { drawLayers(resolveHazardGeom());
  applyFrame(frameSelect.getValue()); applyPlanetFrame(planetFrameSelect.getValue()); });
chkHand.onChange(function ()  { drawLayers(resolveHazardGeom());
  applyFrame(frameSelect.getValue()); applyPlanetFrame(planetFrameSelect.getValue()); });

exportBtn.onClick(function () {
  if (!pending) { return; }
  var tag = 'rasuwa_' + CFG.hazardBasis + '_' + CFG.corridorM + 'm';

  Export.table.toDrive({collection: pending.exposed,
    description: 'exposed_buildings_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});
  Export.table.toDrive({collection: pending.downstream,
    description: 'downstream_reaches_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});
  Export.table.toDrive({
    collection: ee.FeatureCollection([ee.Feature(pending.hazardGeom)]),
    description: 'hazard_footprint_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});
  Export.table.toDrive({collection: pending.disturbVec,
    description: 'observed_disturbance_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'GeoJSON'});

  // summary stats as one-row CSV
  var statRow = ee.Feature(null, {
    event: CFG.eventDate, basis: CFG.hazardBasis, corridor_m: CFG.corridorM,
    n_exposed: pending.exposed.size(),
    people: pop.reduceRegion({reducer: ee.Reducer.sum(),
      geometry: pending.hazardGeom, scale: 100, maxPixels: 1e9,
      bestEffort: true}).get('pop'),
    crop_ha: ee.Number(cropAreaImg.reduceRegion({reducer: ee.Reducer.sum(),
      geometry: pending.hazardGeom, scale: 10, maxPixels: 1e10,
      bestEffort: true, tileScale: 4}).get('carea')).divide(1e4),
    storm_mm: imergTotal.reduceRegion({reducer: ee.Reducer.mean(),
      geometry: basinGeom, scale: imergScale, maxPixels: 1e9,
      bestEffort: true}).get('mm')
  });
  Export.table.toDrive({collection: ee.FeatureCollection([statRow]),
    description: 'summary_stats_' + tag,
    folder: 'nepal_debris_flow', fileFormat: 'CSV'});

  status.setValue('5 tasks queued — open the Tasks tab and press Run on each.');
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
print('Sentinel-1 scenes in change window (pre+post):',
  s1.filterDate(EVENT.advance(-CFG.s1PreDays, 'day'),
                EVENT.advance(CFG.s1PostDays, 'day')).size());

run();
