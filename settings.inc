// Name: settings.inc
// License: Public Domain
// Author: Christopher Bennett
// Version: 1.0, 2025-01-22
// Description: This file contains settings and utility functions for Stellarium scripts

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


/**
 * Sets up the default display settings for the script. The current environment is
 * is saved under the name "RestoreState" before the settings are applied.
 * 
 * The settings applied are:
 * - Clears the current display to "natural" mode.
 *   (az mount, atmosphere, landscape, no lines, labels or markers)
 * - Configures grid lines colors.
 * - Sets font sizes for solar system objects and stars.
 * - Configures star and constellation settings.
 * - Sets observer location and landscape.
 * - Hides the GUI.
 * - Saves the state with the name "defaultState".
 */
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

/**
 * Sets up the display settings for a constellation tour.
 * 
 * The settings applied are:
 * - Clears the current display to "deepspace" mode.
 *   (Eq mount, no atmosphere, landscape etc. No planets, lines, labels or markers)
 * - Sets absolute and relative star scales to 1.0.
 * - Adjusts constellation font size, art intensity, and colors for lines, labels, and boundaries.
 * - Configures constellation line thickness and selection flags.
 * - Enables the display of constellation lines, art, and boundaries.
 * - Saves the current state under the name "constellationState".
 */
function setupConstellations() {
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

/**
 * Cleans up the application state by restoring the previously saved state.
 * It also: 
 * - Manually resets the Solar System and Star font size to their default values (13).
 * - Deletes all markers and labels.
 * - Makes the GUI visible.
 */
function cleanup() {
    restoreState("RestoreState");

    // There is no getter for SolarSystem and StarMgr font size, restore them manually.
    StarMgr.setFontSize(13); // default size is 13
    SolarSystem.setFontSize(13); // default size is 13

    MarkerMgr.deleteAllMarkers();
    LabelMgr.deleteAllLabels();

    core.setGuiVisible(true);
}

/**
 * Sets the time rate in the core system to a very slow rate (0.05) and waits for
 * the user to press the "L" key then restores the original time rate.
 *
 * Note: This function uses a blocking wait loop to effectively pause the script.
 */
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