// Name: utilities.inc
// License: Public Domain
// Author: Christopher Bennett
// Version: 1.0, 2025-01-22
// Description: This file contains utility functions for Stellarium scripts

include("save_state.inc");

// colors
White = core.vec3f(1.0, 1.0, 1.0)                     // #FFFFFF
Gray20 = core.vec3f(0.8, 0.8, 0.8)                    // #CCCCCC
Gray40 = core.vec3f(0.6, 0.6, 0.6)                    // #999999  
Gray60 = core.vec3f(0.4, 0.4, 0.4)                    // #666666
Gray80 = core.vec3f(0.2, 0.2, 0.2)                    // #333333

TrueGreen = core.vec3f(0.00000, 0.56863, 0.00000)     // #009100
BarneyPurple = core.vec3f(0.56863, 0.00000, 0.56863)  // #910091
Cobalt = core.vec3f(0.00000, 0.28627, 0.56863)        // #004991

DarkGreen = core.vec3f(0.12549, 0.37647, 0.12549)     // #206020
DeepViolet = core.vec3f(0.37647, 0.12549, 0.37647)    // #602060
DarkBlue = core.vec3f(0.14118, 0.28627, 0.42745)      // #204160

VeryDarkGreen = core.vec3f(0.17647, 0.28235, 0.17647) // #2D482D
DirtyPurple = core.vec3f(0.28235, 0.17647, 0.28235)   // #482D48
DarkSlate = core.vec3f(0.17647, 0.23137, 0.28235)     // #2D3B48

BloodRed = core.vec3f(0.56863, 0.00000, 0.00392)      // #910001
DarkRed = core.vec3f(0.37647, 0.12549, 0.12941)       // #602021


// Set up the environment for the script
function setup() {
    // Wait... Stellarium needs to start completely.
    core.debug("Setting up environment");

    // reset view to home and save state.
    core.goHome();
    saveState("RestoreState");

    // az mount, atmosphere, landscape, no lines, labels or markers
    core.clear("natural");

    GridLinesMgr.setColorEquatorGrid(DirtyPurple);
    GridLinesMgr.setColorAzimuthalGrid(VeryDarkGreen);

    StelSkyLayerMgr.setFlagShow(false)

    SolarSystem.setFontSize(15); // default size is 13

    StarMgr.setFontSize(15); // default size is 13
    StarMgr.setLabelsAmount(3); // default is 3

    StelSkyDrawer.setAbsoluteStarScale(1.0);
    StelSkyDrawer.setRelativeStarScale(1.0);
    ZodiacalLight.setFlagShow(false);

    StelMovementMgr.setAutoMoveDuration(3);

    ConstellationMgr.setFontSize(18); // default size is 15
    ConstellationMgr.setArtIntensity(0.2);
    ConstellationMgr.setLinesColor(Gray40);
    ConstellationMgr.setLabelsColor(Gray40);
    ConstellationMgr.setBoundariesColor(Gray80);
    ConstellationMgr.setConstellationLineThickness(1);
    ConstellationMgr.setFlagConstellationPick(true);
    ConstellationMgr.setFlagIsolateSelected(true);
    ConstellationMgr.deselectConstellations();

    core.wait(2);
    // core.setObserverLocation("Bear Branch", "Earth");
    core.setObserverLocation(-76.986769, 39.647307, 210, 2, "BBNC", "Earth");
    LandscapeMgr.setCurrentLandscapeName("Bear Branch");

    core.setGuiVisible(false);

    saveState("defaultState");
}

// Set up the environment for constellations tour
function setupConstellations() {
    // Eq mount, no atmosphere, landscape etc. No planets, lines, labels or markers.
    core.clear("deepspace");

    StelSkyDrawer.setAbsoluteStarScale(1.0);
    StelSkyDrawer.setRelativeStarScale(1.0);

    ConstellationMgr.setFontSize(18); // default size is 15
    ConstellationMgr.setArtIntensity(0.2);
    ConstellationMgr.setLinesColor(Gray40);
    ConstellationMgr.setLabelsColor(Gray40);
    ConstellationMgr.setBoundariesColor(Gray80);
    ConstellationMgr.setConstellationLineThickness(1);
    ConstellationMgr.setFlagConstellationPick(true);
    ConstellationMgr.setFlagIsolateSelected(true);
    ConstellationMgr.deselectConstellations();

    ConstellationMgr.setFlagLines(true);
    ConstellationMgr.setFlagArt(true);
    ConstellationMgr.setFlagBoundaries(true);

    saveState("constellationState");
}

// Try to undo script settings that will mess up Stellarium's expected operation
function cleanup() {
    // reset view to home.
    restoreState("RestoreState");

    // There is no getter for SolarSystem and StarMgr font size, restore them manually.
    StarMgr.setFontSize(13); // default size is 13
    SolarSystem.setFontSize(13); // default size is 13

    MarkerMgr.deleteAllMarkers();
    LabelMgr.deleteAllLabels();

    core.setGuiVisible(true);
}

// Pause the script until key press "L".
function pauseKey() {
    var timerate = core.getTimeRate();
    core.setTimeRate(0.05);
    i = 0;
    do {
        core.wait(1);
        i++;
    }
    while (core.getTimeRate() < 0.1);
    core.setTimeRate(timerate);
}