/**** Conditioning factors for the Dhaka / Chattogram subsidence XAI paper ****
 *
 * Run in the Earth Engine Code Editor (code.earthengine.google.com).
 * Every layer is aggregated onto the exact 100 m lon/lat lattice of the
 * three-city decomposition (same cells as <city>_factors.csv), and exported as
 * CSV with lattice row/col so the merge is exact. Ten table exports and three
 * GeoTIFFs per city (26 tasks in all) keep each task small; start them all from the
 * Tasks tab.
 *
 *   A_ob1/2   Open Buildings 2.5D Temporal, read at 4 m (its effective resolution), two tasks
 *   A_urban   GHSL built surface/volume/height and urbanisation epoch, WorldPop, VIIRS
 *   B_gsw     JRC GSW history (fill year, water lost during the InSAR window)
 *   B_ls      Landsat 5 water 1990/2005 and Landsat 8 water 2022 (same method,
 *             with observation counts)
 *   B_dw      Dynamic World water/built at both ends of the InSAR window
 *   C_s2      Sentinel-2 dry-season indices
 *   C_lst     Landsat 8/9 dry-season land surface temperature
 *   C_dw      Dynamic World class probabilities 2017-2023
 *   D_soil    OpenLandMap clay/sand/bulk density/organic carbon/texture, MERIT Hydro
 *             (CHIRPS and SMAP are not exported: at 5-9 km they are dropped by the
 *             1 km resolution rule anyway, and the SMAP collection version changes)
 *   E_dem     GeoTIFFs at 30 m UTM 46N: Copernicus GLO-30 and FABDEM (bilinear), and
 *             the GSW permanent-water mask for the distance to water
 *
 * Output folder on Drive: XAI_factors.  Send all CSVs and the three tifs back
 * (<c>_glo30.tif, <c>_fabdem.tif, <c>_gsw_water.tif).
 * Set USE_COMMUNITY = false if a community-catalogue asset (FABDEM) fails to load.
 * v2 (18 Oct 2026): urban epoch without the 2030 fill value, water loss from one
 * method, GSW loss inside the InSAR window, height change of existing buildings only,
 * bilinear DEMs, Open Buildings split out, Landsat observation counts, fuller cloud masks.
 * v3 (19 Oct 2026): 6 km margin set in degrees (no buffer), empty-collection guard for Landsat, point
 * geometry in every CSV so the merge can check row/col, OB_SPLIT fallback, larger
 * Chattogram lattice (whole city), GSW window loss described as 2017-2019.
 *
 * v4 (19 Oct 2026): "User memory limit exceeded" in dha_A_ob and dha_C_surface. Open
 * Buildings is now read at 4 m (its effective resolution; the 0.5 m grid is an upsampled
 * product) and exported in two tasks; C is split into C_s2 and C_dw, the Dynamic World
 * label mode is dropped (the class probabilities carry it), and every table export runs
 * with tileScale 16. If an Open Buildings task still runs out of memory, set
 * OB_SPLIT = true: each band then becomes its own task.
 * v5 (19 Oct 2026): B_fill and C_s2 failed within a minute. B is split into B_gsw, B_ls and
 * B_dw, and C_s2 into C_s2 and C_lst, so a failure points at one source; Sentinel-2 uses
 * four dry seasons, scenes under 50 % cloud, five bands and a 20 m composite; Landsat skips
 * scenes over 70 % cloud; the empty-epoch guard no longer merges a placeholder into the
 * collection but swaps it in only when an epoch has no scenes.
 */

var GRIDS = {
 "dha": {
  "a": 90.26159662822207,
  "b": 0.000982773746135546,
  "c": 23.552355546646297,
  "d": 0.0008992940037785201,
  "r0": 106,
  "r1": 437,
  "c0": 55,
  "c1": 328,
  "nrows": 332,
  "ncols": 274,
  "west": 90.31515779738645,
  "east": 90.5844378038276,
  "south": 23.64723106404493,
  "north": 23.945796673299398,
  "crs": "EPSG:4326",
  "crsTransform": [
   0.000982773746135546,
   0,
   90.31515779738645,
   0,
   -0.0008992940037785201,
   23.945796673299398
  ]
 },
 "ctg": {
  "a": 91.7391702201778,
  "b": 0.0009722341308043657,
  "c": 22.186885097352683,
  "d": 0.000899277134176499,
  "r0": 10,
  "r1": 314,
  "c0": 2,
  "c1": 176,
  "nrows": 305,
  "ncols": 175,
  "west": 91.740628571374,
  "east": 91.91076954426477,
  "south": 22.19542823012736,
  "north": 22.469707756051193,
  "crs": "EPSG:4326",
  "crsTransform": [
   0.0009722341308043657,
   0,
   91.740628571374,
   0,
   -0.000899277134176499,
   22.469707756051193
  ]
 }
};            // filled from data/<c>_grid.json
var FOLDER = 'XAI_factors';
var USE_COMMUNITY = true;
var OB_SPLIT = false;
// which exports this run creates (true = create the task)
var EXPORT = {"A_ob": false, "A_urban": false, "B_gsw": true, "B_ls": true, "B_dw": true, "C_s2": true, "C_lst": true, "C_dw": false, "D_soil": false, "E_dem": false};
var T0 = '2017-10-01', T1 = '2023-10-01';   // InSAR window (11 Oct 2017 - 20 Sep 2023)

function dryMonths() {             // Nov-Mar: avoids monsoon flooding and clouds
  return ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));
}

function run(city) {
  var g = GRIDS[city];
  var proj = ee.Projection(g.crs, g.crsTransform);
  var region = ee.Geometry.Rectangle([g.west, g.south, g.east, g.north], null, false);
  // 6 km margin added directly in degrees (no buffer/bounds geometry operations, which need
  // error margins); it covers the 5 km work grids of s32a/s32c even at the corners
  var lat0 = (g.south + g.north) / 2;
  var mdy = 6000 / 111320, mdx = 6000 / (111320 * Math.cos(lat0 * Math.PI / 180));
  var work = ee.Geometry.Rectangle([g.west - mdx, g.south - mdy, g.east + mdx, g.north + mdy], null, false);

  // fine image -> 100 m lattice; categorical layers use mode
  function toGrid(img, reducer, maxPix) {
    return img.reduceResolution({reducer: reducer || ee.Reducer.mean(), maxPixels: maxPix || 4096})
              .reproject(proj);
  }
  function via4(img) {   // Open Buildings: read at 4 m (nearest), then 4 m -> 100 m lattice mean
    return toGrid(img.reproject(ee.Projection('EPSG:32646').atScale(4)), ee.Reducer.mean(), 1024);
  }
  // lattice indices from pixel-centre lon/lat: col = round((lon - a)/b), row = round((lat - c)/d)
  var ll = ee.Image.pixelLonLat();
  var rowcol = ll.select('latitude').subtract(g.c).divide(g.d).round().rename('row')
      .addBands(ll.select('longitude').subtract(g.a).divide(g.b).round().rename('col')).toInt();

  // geometries: true writes each sampled pixel centre to .geo, so s32h can recompute row/col
  // from the coordinates and check them against the row/col bands
  function exportTable(img, name, tileScale) {
    var key = name.indexOf('A_ob') === 0 ? 'A_ob' : name;
    if (!EXPORT[key]) { return; }
    var t = img.addBands(rowcol).reproject(proj).sample({
      region: region, projection: proj, dropNulls: false, geometries: true, tileScale: tileScale || 16});
    Export.table.toDrive({collection: t, description: city + '_' + name, folder: FOLDER,
                          fileNamePrefix: city + '_' + name, fileFormat: 'CSV'});
  }

  // ------------------------------------------------------------------ A urban load
  var ob = ee.ImageCollection('GOOGLE/Research/open-buildings-temporal/v1').filterBounds(work);
  function obYear(y) {
    var col = ob.filterDate(y + '-01-01', (y + 1) + '-01-01');
    return col.mosaic().setDefaultProjection(col.first().projection());
  }
  var ob16 = obYear(2016), ob23 = obYear(2023);
  var pres16 = ob16.select('building_presence').gt(0.5);
  var pres23 = ob23.select('building_presence').gt(0.5);
  var area = ee.Image.pixelArea();
  var AOB1 = [
    via4(pres16.rename('ob_bfrac_2016')),
    via4(pres23.rename('ob_bfrac_2023')),
    via4(pres23.and(pres16.not()).rename('ob_new_frac'))                    // new footprint 2016-2023
  ];
  var AOB2 = [
    via4(ob23.select('building_height').updateMask(pres23).rename('ob_height_2023')),
    via4(ob23.select('building_height').multiply(pres23).rename('ob_volume_2023')),   // m3 per m2
    // height change of buildings present in both years; new buildings are in ob_new_frac
    via4(ob23.select('building_height').subtract(ob16.select('building_height'))
             .updateMask(pres23.and(pres16)).rename('ob_dheight_existing')),
    via4(ob23.select('building_fractional_count').divide(area).multiply(1e4).rename('ob_bld_per_ha'))
  ];
  if (OB_SPLIT) {
    AOB1.concat(AOB2).forEach(function (im, i) { exportTable(ee.Image(im).toFloat(), 'A_ob_b' + i); });
  } else {
    exportTable(ee.Image.cat(AOB1).toFloat(), 'A_ob1');
    exportTable(ee.Image.cat(AOB2).toFloat(), 'A_ob2');
  }
  var epochs = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020];
  var bs = function (y) { return ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/' + y).select('built_surface'); };
  // first epoch in which a 100 m GHSL pixel is >= 20 % built (2000 m2); never-built pixels stay
  // masked (no 2030 fill value, which the mean would mix with real years); min over the pixels
  // that touch a lattice cell; ghs_urban_by2020 = share of the cell's pixels urban by 2020
  var urbEpoch = ee.ImageCollection(epochs.map(function (y) {
    return ee.Image.constant(y).updateMask(bs(y).gte(2000)).toFloat();
  })).min().rename('ghs_urban_epoch');
  var ghsProj = bs(2020).projection();
  var A = ee.Image.cat([
    toGrid(urbEpoch.setDefaultProjection(ghsProj), ee.Reducer.min()),
    toGrid(bs(2020).gte(2000).rename('ghs_urban_by2020')),
    toGrid(bs(1990).rename('ghs_bs_1990')), toGrid(bs(2020).rename('ghs_bs_2020')),
    toGrid(bs(2020).subtract(bs(2000)).rename('ghs_bs_d2000_2020')),
    toGrid(ee.Image('JRC/GHSL/P2023A/GHS_BUILT_V/2020').select('built_volume_total').rename('ghs_bv_2020')),
    toGrid(ee.Image('JRC/GHSL/P2023A/GHS_BUILT_H/2018').select('built_height').rename('ghs_bh_2018'))
  ]);
  var wp = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('country', 'BGD'))
             .filter(ee.Filter.eq('year', 2020)).first();
  A = A.addBands(toGrid(wp.select('population').divide(ee.Image.pixelArea()).multiply(1e6)
                         .setDefaultProjection(wp.projection()).rename('pop_density_2020')));
  var viirs = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG').filterDate(T0, T1).select('avg_rad');
  A = A.addBands(toGrid(viirs.median().setDefaultProjection(viirs.first().projection()).rename('viirs_rad')));
  exportTable(A.toFloat(), 'A_urban');

  // ------------------------------------------------------------------ B fill / reclamation
  var yh = ee.ImageCollection('JRC/GSW1_4/YearlyHistory').filterBounds(work);
  var gswProj = yh.first().projection();
  var wetYear = yh.map(function (im) {
    var y = ee.Number.parse(ee.Date(im.get('system:time_start')).format('YYYY'));
    return ee.Image.constant(y).updateMask(im.select('waterClass').gte(2)).toFloat().rename('y');
  });
  var everWet = yh.map(function (im) { return im.select('waterClass').gte(2); }).max();
  var recentWet = yh.filterDate('2019-01-01', '2022-01-01')
                    .map(function (im) { return im.select('waterClass').gte(2); }).max();
  var lastWet = wetYear.max();
  var lost = everWet.and(recentWet.not());                   // water once, dry in 2019-2021 = filled/drained
  var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
  var B = ee.Image.cat([          // B_gsw: water history from JRC GSW
    toGrid(gsw.select('occurrence').unmask(0).rename('gsw_occurrence')),
    toGrid(everWet.unmask(0).rename('gsw_ever_water').setDefaultProjection(gswProj)),
    toGrid(lost.unmask(0).rename('gsw_lost_share').setDefaultProjection(gswProj)),
    toGrid(lastWet.add(1).updateMask(lost).rename('gsw_fill_year').setDefaultProjection(gswProj)),
    toGrid(lost.and(lastWet.gte(2010)).unmask(0).rename('gsw_lost_after2010').setDefaultProjection(gswProj)),
    // last wet year 2016-2018 = water lost in 2017-2019, inside the InSAR window ('lost' needs the
    // pixel dry through 2019-2021 and GSW ends in 2021, so later losses cannot be seen).
    // 'lost' means drained OR filled (embanked floodplain counts too): read with dw_built_change / ob_new_frac
    toGrid(lost.and(lastWet.gte(2016)).unmask(0).rename('gsw_lost_window').setDefaultProjection(gswProj))
  ]);
  exportTable(B.toFloat(), 'B_gsw');
  // Landsat dry-season MNDWI (green, SWIR1), one method for every epoch, so water loss is
  // water_share_1990 - water_share_2022. Clear = not fill/dilated cloud/cirrus/cloud/shadow.
  function clearLandsat(qa) {
    return qa.bitwiseAnd(1 << 0).or(qa.bitwiseAnd(1 << 1)).or(qa.bitwiseAnd(1 << 2))
             .or(qa.bitwiseAnd(1 << 3)).or(qa.bitwiseAnd(1 << 4)).eq(0);
  }
  // guard: if an epoch has no scenes (the late-1980s Landsat 5 archive over Bangladesh is thin)
  // use a fully masked image instead of failing; it never adds a value
  function orEmpty(col, name) {
    var blank = ee.Image.constant(0).toFloat().rename(name).updateMask(0);
    return ee.ImageCollection(ee.Algorithms.If(col.size().gt(0), col, ee.ImageCollection([blank])));
  }
  function lsMndwi(id, green, swir, y0, y1) {
    var c = orEmpty(ee.ImageCollection(id).filterBounds(work)
      .filterDate(y0 + '-01-01', y1 + '-12-31').filter(dryMonths())
      .filter(ee.Filter.lt('CLOUD_COVER', 70))                 // skip near-cloudy scenes: less to hold
      .map(function (im) {
        var sr = im.select([green, swir]).multiply(0.0000275).add(-0.2);
        return sr.normalizedDifference([green, swir]).rename('mndwi').updateMask(clearLandsat(im.select('QA_PIXEL')));
      }), 'mndwi');
    var p30 = ee.Projection('EPSG:32646').atScale(30);
    return {m: c.median().setDefaultProjection(p30), n: c.count().setDefaultProjection(p30)};
  }
  var l90 = lsMndwi('LANDSAT/LT05/C02/T1_L2', 'SR_B2', 'SR_B5', 1988, 1992);
  var l05 = lsMndwi('LANDSAT/LT05/C02/T1_L2', 'SR_B2', 'SR_B5', 2003, 2007);
  var l22 = lsMndwi('LANDSAT/LC08/C02/T1_L2', 'SR_B3', 'SR_B6', 2021, 2023);
  var BLS = ee.Image.cat([        // B_ls: Landsat water 1990 / 2005 / 2022
    toGrid(l90.m.rename('mndwi_1990')), toGrid(l90.m.gt(0).rename('water_share_1990')),
    toGrid(l90.n.rename('l5_nobs_1990')),
    toGrid(l05.m.rename('mndwi_2005')), toGrid(l05.m.gt(0).rename('water_share_2005')),
    toGrid(l05.n.rename('l5_nobs_2005')),
    toGrid(l22.m.rename('mndwi_2022')), toGrid(l22.m.gt(0).rename('water_share_2022')),
    toGrid(l22.n.rename('l8_nobs_2022'))
  ]);
  exportTable(BLS.toFloat(), 'B_ls');
  var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filterBounds(work);
  var dwProj = dw.first().projection();
  function dwMean(band, a, b) { return dw.filterDate(a, b).select(band).mean().setDefaultProjection(dwProj); }
  var BDW = ee.Image.cat([        // B_dw: Dynamic World built/water at both ends of the window
    toGrid(dwMean('water', '2017-10-01', '2018-10-01').rename('dw_water_2018')),
    toGrid(dwMean('built', '2017-10-01', '2018-10-01').rename('dw_built_2018')),
    toGrid(dwMean('built', '2022-10-01', '2023-10-01').rename('dw_built_2023'))
  ]);
  exportTable(BDW.toFloat(), 'B_dw');

  // ------------------------------------------------------------------ C surface (S2, DW, LST)
  var cs = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
  // four dry seasons inside the InSAR window, scenes under 50 % cloud, only the bands needed,
  // composited at 20 m (the SWIR bands are 20 m anyway; 25 pixels per 100 m cell is plenty)
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(work)
    .filterDate('2019-11-01', '2023-04-01').filter(dryMonths())
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 50))
    .linkCollection(cs, ['cs_cdf'])
    .map(function (im) {
      return im.select(['B2', 'B3', 'B4', 'B8', 'B11']).updateMask(im.select('cs_cdf').gte(0.6)).divide(10000);
    });
  var s2m = s2.median().setDefaultProjection(ee.Projection('EPSG:32646').atScale(20));
  var ndvi = s2m.normalizedDifference(['B8', 'B4']);
  var C = ee.Image.cat([          // C_s2: Sentinel-2 indices
    toGrid(ndvi.rename('s2_ndvi')),
    toGrid(s2m.normalizedDifference(['B11', 'B8']).rename('s2_ndbi')),
    toGrid(s2m.normalizedDifference(['B3', 'B11']).rename('s2_mndwi')),
    toGrid(s2m.expression('((b("B11")+b("B4"))-(b("B8")+b("B2")))/((b("B11")+b("B4"))+(b("B8")+b("B2")))').rename('s2_bsi')),
    toGrid(ndvi.rename('s2_ndvi_sd'), ee.Reducer.stdDev())
  ]);
  // C_dw: Dynamic World class probabilities over the InSAR window (one band list, one mean)
  var dwb = ['built', 'trees', 'grass', 'crops', 'bare', 'water', 'flooded_vegetation'];
  var CDW = toGrid(dw.filterDate(T0, T1).select(dwb).mean().setDefaultProjection(dwProj)
                     .rename(dwb.map(function (b) { return 'dw_' + b; })));
  exportTable(CDW.toFloat(), 'C_dw');
  exportTable(C.toFloat(), 'C_s2');
  // C_lst: Landsat 8/9 dry-season land surface temperature (context only)
  var lst = orEmpty(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(work).filterDate('2018-01-01', '2024-01-01').filter(dryMonths())
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(function (im) {
      return im.select('ST_B10').multiply(0.00341802).add(149.0 - 273.15).rename('lst')
               .updateMask(clearLandsat(im.select('QA_PIXEL')));
    }), 'lst').median().setDefaultProjection(ee.Projection('EPSG:32646').atScale(30));
  exportTable(toGrid(lst.rename('lst_dry')).toFloat(), 'C_lst');

  // ------------------------------------------------------------------ D soil, hydrology, climate
  var olm = function (id) { return ee.Image('OpenLandMap/SOL/' + id); };
  // depth-weighted 0-200 cm mean from the six standard depths (trapezoid rule)
  function dmean(img) {
    var d = [0, 10, 30, 60, 100, 200], b = ['b0', 'b10', 'b30', 'b60', 'b100', 'b200'];
    var s = ee.Image(0);
    for (var i = 0; i < 5; i++) {
      s = s.add(img.select(b[i]).add(img.select(b[i + 1])).divide(2).multiply(d[i + 1] - d[i]));
    }
    return s.divide(200).setDefaultProjection(img.projection());
  }
  var clay = olm('SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'), sand = olm('SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02');
  var bd = olm('SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02'), oc = olm('SOL_ORGANIC-CARBON_USDA-6A1C_M/v02');
  var tex = olm('SOL_TEXTURE-CLASS_USDA-TT_M/v02');
  var mh = ee.Image('MERIT/Hydro/v1_0_1');
  var D = ee.Image.cat([
    toGrid(dmean(clay).rename('olm_clay')), toGrid(clay.select('b0').rename('olm_clay_top')),
    toGrid(dmean(sand).rename('olm_sand')),
    toGrid(dmean(bd).multiply(10).rename('olm_bulkdens')),          // kg/m3
    toGrid(dmean(oc).multiply(5).rename('olm_orgc')),               // g/kg
    toGrid(tex.select('b30').rename('olm_texture_30cm'), ee.Reducer.mode()),
    toGrid(mh.select('hnd').rename('merit_hand')),
    toGrid(mh.select('upa').log10().rename('merit_log_upa')),
    toGrid(mh.select('elv').rename('merit_elev'))
  ]);
  exportTable(D.toFloat(), 'D_soil');

  // ------------------------------------------------------------------ E DEMs (30 m UTM GeoTIFF)
  if (!EXPORT.E_dem) { Map.addLayer(region, {color: 'red'}, city + ' lattice'); return; }
  // bilinear: the default nearest-neighbour step would be resampled a second time in s32a
  var glo = ee.ImageCollection('COPERNICUS/DEM/GLO30').filterBounds(work).select('DEM');
  Export.image.toDrive({image: glo.mosaic().setDefaultProjection(glo.first().projection()).resample('bilinear').toFloat(),
    description: city + '_E_glo30', folder: FOLDER, fileNamePrefix: city + '_glo30', region: work,
    crs: 'EPSG:32646', scale: 30, maxPixels: 1e9});
  if (USE_COMMUNITY) {
    var fab = ee.ImageCollection('projects/sat-io/open-datasets/FABDEM').filterBounds(work);
    Export.image.toDrive({image: fab.mosaic().setDefaultProjection(fab.first().projection()).resample('bilinear').toFloat(),
      description: city + '_E_fabdem', folder: FOLDER, fileNamePrefix: city + '_fabdem', region: work,
      crs: 'EPSG:32646', scale: 30, maxPixels: 1e9});
  }
  // permanent surface water (GSW occurrence >= 50 %), full coverage, for dist_water in s32c
  Export.image.toDrive({image: gsw.select('occurrence').unmask(0).gte(50).toByte(),
    description: city + '_E_gsw_water', folder: FOLDER, fileNamePrefix: city + '_gsw_water', region: work,
    crs: 'EPSG:32646', scale: 30, maxPixels: 1e9});
  Map.addLayer(region, {color: 'red'}, city + ' lattice');
}

run('dha');
run('ctg');
Map.centerObject(ee.Geometry.Point([91.0, 23.1]), 8);
