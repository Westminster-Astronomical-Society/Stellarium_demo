// Name: solar_calc.js
// License: Public Domain
// Author: Chris Bennett
// Version: 1.0
// Description: Script for Solar calculations. Mostly stolen from
// http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html


/* Solar position calculation functions */
/*************************************************************/

/**
 * Converts an angle from radians to degrees.
 *
 * @param {number} angleRad - The angle in radians.
 * @returns {number} The angle in degrees.
 */
function radToDeg(angleRad) {
    return (180.0 * angleRad / Math.PI);
}

/**
 * Converts an angle from degrees to radians.
 *
 * @param {number} angleDeg - The angle in degrees.
 * @returns {number} The angle in radians.
 */
function degToRad(angleDeg) {
    return (Math.PI * angleDeg / 180.0);
}

/**
 * Calculates the Julian Century from the Julian Date.
 *
 * @param {number} jd - The Julian Date.
 * @returns {number} The Julian Century.
 */
function calcTimeJulianCent(jd) {
    var T = (jd - 2451545.0) / 36525.0
    return T
}

/**
 * Calculates the Julian Date (JD) from Julian Centuries (t).
 *
 * @param {number} t - The number of Julian Centuries since J2000.0.
 * @returns {number} The Julian Date corresponding to the given Julian Centuries.
 */
function calcJDFromJulianCent(t) {
    var JD = t * 36525.0 + 2451545.0
    return JD
}

/**
 * Determines if a given year is a leap year.
 *
 * A leap year is exactly divisible by 4 except for end-of-century years, which must be divisible by 400.
 * This means that the year 2000 was a leap year, although 1900 was not.
 *
 * @param {number} yr - The year to be checked.
 * @returns {boolean} - Returns true if the year is a leap year, otherwise false.
 */
function isLeapYear(yr) {
    return ((yr % 4 == 0 && yr % 100 != 0) || yr % 400 == 0);
}

/**
 * Converts a Julian Date (JD) to a calendar date.
 *
 * @param {number} jd - The Julian Date to convert.
 * @returns {Object} An object containing the year, month, and day corresponding to the given Julian Date.
 * @returns {number} return.year - The year of the calendar date.
 * @returns {number} return.month - The month of the calendar date (1-12).
 * @returns {number} return.day - The day of the calendar date.
 */
function calcDateFromJD(jd) {
    var z = Math.floor(jd + 0.5);
    var f = (jd + 0.5) - z;
    if (z < 2299161) {
        var A = z;
    } else {
        alpha = Math.floor((z - 1867216.25) / 36524.25);
        var A = z + 1 + alpha - Math.floor(alpha / 4);
    }
    var B = A + 1524;
    var C = Math.floor((B - 122.1) / 365.25);
    var D = Math.floor(365.25 * C);
    var E = Math.floor((B - D) / 30.6001);
    var day = B - D - Math.floor(30.6001 * E) + f;
    var month = (E < 14) ? E - 1 : E - 13;
    var year = (month > 2) ? C - 4716 : C - 4715;

    return { "year": year, "month": month, "day": day }
}

/**
 * Calculates the day of the year (DOY) from a given Julian Date (JD).
 *
 * @param {number} jd - The Julian Date.
 * @returns {number} The day of the year (1-366).
 */
function calcDoyFromJD(jd) {
    var date = calcDateFromJD(jd)

    var k = (isLeapYear(date.year) ? 1 : 2);
    var doy = Math.floor((275 * date.month) / 9) - k * Math.floor((date.month + 9) / 12) + date.day - 30;

    return doy;
}

/**
 * Calculates the Julian Date (JD) for a given Gregorian calendar date.
 *
 * @param {number} year - The year of the date.
 * @param {number} month - The month of the date (1-12).
 * @param {number} day - The day of the date (1-31).
 * @returns {number} The Julian Date corresponding to the given date.
 */
function getJD(year, month, day) {
    if (month <= 2) {
        year -= 1
        month += 12
    }
    var A = Math.floor(year / 100)
    var B = 2 - A + Math.floor(A / 4)
    var JD = Math.floor(365.25 * (year + 4716)) + Math.floor(30.6001 * (month + 1)) + day + B - 1524.5
    return JD
}

/**
 * Calculates the geometric mean longitude of the sun.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The geometric mean longitude of the sun in degrees.
 */
function calcGeomMeanLongSun(t) {
    var L0 = 280.46646 + t * (36000.76983 + t * (0.0003032))
    while (L0 > 360.0) {
        L0 -= 360.0
    }
    while (L0 < 0.0) {
        L0 += 360.0
    }
    return L0
}

/**
 * Calculate the geometric mean anomaly of the Sun.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The geometric mean anomaly of the Sun in degrees.
 */
function calcGeomMeanAnomalySun(t) {
    var M = 357.52911 + t * (35999.05029 - 0.0001537 * t);
    return M;
}

/**
 * Calculates the eccentricity of Earth's orbit.
 *
 * The eccentricity of Earth's orbit changes over time due to gravitational interactions with other bodies in the solar system.
 *
 * @param {number} t - Time in Julian centuries since J2000.0.
 * @returns {number} The eccentricity of Earth's orbit (unitless).
 */
function calcEccentricityEarthOrbit(t) {
    var e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t);
    return e;
}

/**
 * Calculate the Sun's equation of center.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The Sun's equation of center in degrees.
 */
function calcSunEqOfCenter(t) {
    var m = calcGeomMeanAnomalySun(t);
    var mrad = degToRad(m);
    var sinm = Math.sin(mrad);
    var sin2m = Math.sin(mrad + mrad);
    var sin3m = Math.sin(mrad + mrad + mrad);
    var C = sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289;
    return C;
}

/**
 * Calculates the true longitude of the sun.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The true longitude of the sun in degrees.
 */
function calcSunTrueLong(t) {
    var l0 = calcGeomMeanLongSun(t);
    var c = calcSunEqOfCenter(t);
    var O = l0 + c;
    return O;
}

/**
 * Calculates the Sun's true anomaly.
 *
 * @param {number} t - The time in Julian centuries since J2000.0.
 * @returns {number} The Sun's true anomaly in degrees.
 */
function calcSunTrueAnomaly(t) {
    var m = calcGeomMeanAnomalySun(t);
    var c = calcSunEqOfCenter(t);
    var v = m + c;
    return v;
}

/**
 * Calculates the Sun's radius vector (distance from the Earth to the Sun) in astronomical units (AU).
 *
 * @param {number} t - The time in Julian centuries since J2000.0.
 * @returns {number} The Sun's radius vector in astronomical units (AU).
 */
function calcSunRadVector(t) {
    var v = calcSunTrueAnomaly(t);
    var e = calcEccentricityEarthOrbit(t);
    var R = (1.000001018 * (1 - e * e)) / (1 + e * Math.cos(degToRad(v)));
    return R;
}

/**
 * Calculates the apparent longitude of the Sun.
 *
 * @param {number} t - Julian century since J2000.0.
 * @returns {number} The apparent longitude of the Sun in degrees.
 */
function calcSunApparentLong(t) {
    var o = calcSunTrueLong(t);
    var omega = 125.04 - 1934.136 * t;
    var lambda = o - 0.00569 - 0.00478 * Math.sin(degToRad(omega));
    return lambda;
}

/**
 * Calculates the mean obliquity of the ecliptic.
 *
 * The obliquity of the ecliptic is the angle between the plane of the Earth's orbit and the plane of the Earth's equator.
 * This function uses a formula to approximate the mean obliquity of the ecliptic at a given time.
 *
 * @param {number} t - Julian centuries since J2000.0.
 * @returns {number} The mean obliquity of the ecliptic in degrees.
 */
function calcMeanObliquityOfEcliptic(t) {
    var seconds = 21.448 - t * (46.8150 + t * (0.00059 - t * (0.001813)));
    var e0 = 23.0 + (26.0 + (seconds / 60.0)) / 60.0;
    return e0;
}

/**
 * Calculates the obliquity correction for a given time.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The corrected obliquity of the ecliptic in degrees.
 */
function calcObliquityCorrection(t) {
    var e0 = calcMeanObliquityOfEcliptic(t);
    var omega = 125.04 - 1934.136 * t;
    var e = e0 + 0.00256 * Math.cos(degToRad(omega));
    return e;
}

/**
 * Calculates the Right Ascension (RA) of the Sun.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The Right Ascension of the Sun in degrees.
 */
function calcSunRA(t) {
    var e = calcObliquityCorrection(t);
    var lambda = calcSunApparentLong(t);
    var tananum = (Math.cos(degToRad(e)) * Math.sin(degToRad(lambda)));
    var tanadenom = (Math.cos(degToRad(lambda)));
    var alpha = radToDeg(Math.atan2(tananum, tanadenom));
    return alpha;
}

/**
 * Calculates the declination of the Sun.
 *
 * @param {number} t - The number of Julian centuries since J2000.0.
 * @returns {number} The declination of the Sun in degrees.
 */
function calcSunDec(t) {
    var e = calcObliquityCorrection(t);
    var lambda = calcSunApparentLong(t);
    var sint = Math.sin(degToRad(e)) * Math.sin(degToRad(lambda));
    var theta = radToDeg(Math.asin(sint));
    return theta;
}

/**
 * Calculates the equation of time for a given Julian century.
 * The equation of time is the difference between apparent solar time and mean solar time.
 *
 * @param {number} t - Julian centuries since J2000.0.
 * @returns {number} - The equation of time in minutes of time.
 */
function calcEquationOfTime(t) {
    var epsilon = calcObliquityCorrection(t);
    var l0 = calcGeomMeanLongSun(t);
    var e = calcEccentricityEarthOrbit(t);
    var m = calcGeomMeanAnomalySun(t);

    var y = Math.tan(degToRad(epsilon) / 2.0);
    y *= y;

    var sin2l0 = Math.sin(2.0 * degToRad(l0));
    var sinm = Math.sin(degToRad(m));
    var cos2l0 = Math.cos(2.0 * degToRad(l0));
    var sin4l0 = Math.sin(4.0 * degToRad(l0));
    var sin2m = Math.sin(2.0 * degToRad(m));

    var Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;
    return radToDeg(Etime) * 4.0;
}

/**
 * Calculates the hour angle of the sun at sunset or sunrise.
 * (for sunset, use -HA)
 *
 * @param {number} lat - The latitude in degrees.
 * @param {number} solarDec - The solar declination in degrees.
 * @param {number} horizon - The zenith angle of the horizon in degrees.
 * @returns {number} The hour angle of the sun in radians.
 */
function calcHourAngleSun(lat, solarDec, horizon) {
    var latRad = degToRad(lat);
    var sdRad = degToRad(solarDec);
    var HAarg = (Math.cos(degToRad(horizon)) / (Math.cos(latRad) * Math.cos(sdRad)) - Math.tan(latRad) * Math.tan(sdRad));
    var HA = Math.acos(HAarg);
    return HA;
}

/**
 * Checks if the given input is a valid number.
 * The input can be a positive or negative number and can contain one decimal point.
 *
 * @param {string|number} inputVal - The value to be checked.
 * @returns {boolean} - Returns true if the input is a valid number, otherwise false.
 */
function isNumber(inputVal) {
    var oneDecimal = false;
    var inputStr = "" + inputVal;
    for (var i = 0; i < inputStr.length; i++) {
        var oneChar = inputStr.charAt(i);
        if (i == 0 && (oneChar == "-" || oneChar == "+")) {
            continue;
        }
        if (oneChar == "." && !oneDecimal) {
            oneDecimal = true;
            continue;
        }
        if (oneChar < "0" || oneChar > "9") {
            return false;
        }
    }
    return true;
}

/**
 * Calculates the atmospheric refraction correction for a given elevation angle.
 *
 * @param {number} elev - The elevation angle in degrees.
 * @returns {number} The refraction correction in degrees.
 */
function calcRefraction(elev) {

    if (elev > 85.0) {
        var correction = 0.0;
    } else {
        var te = Math.tan(degToRad(elev));
        if (elev > 5.0) {
            var correction = 58.1 / te - 0.07 / (te * te * te) + 0.000086 / (te * te * te * te * te);
        } else if (elev > -0.575) {
            var correction = 1735.0 + elev * (-518.2 + elev * (103.4 + elev * (-12.79 + elev * 0.711)));
        } else {
            var correction = -20.774 / te;
        }
        correction = correction / 3600.0;
    }

    return correction
}

/**
 * Calculates the azimuth and elevation of the sun for a given time and location.
 *
 * @param {number} T - Julian century since J2000.0.
 * @param {number} localtime - Local time in minutes from midnight.
 * @param {number} latitude - Latitude of the observer in degrees.
 * @param {number} longitude - Longitude of the observer in degrees.
 * @param {number} zone - Time zone offset from UTC in hours.
 * @returns {Object} An object containing the azimuth and elevation of the sun.
 * @returns {number} azimuth - The azimuth angle of the sun in degrees.
 * @returns {number} elevation - The elevation angle of the sun in degrees.
 */
function calcAzEl(T, localtime, latitude, longitude, zone) {

    var eqTime = calcEquationOfTime(T)
    var theta = calcSunDec(T)

    var solarTimeFix = eqTime + 4.0 * longitude - 60.0 * zone
    var earthRadVec = calcSunRadVector(T)
    var trueSolarTime = localtime + solarTimeFix
    while (trueSolarTime > 1440) {
        trueSolarTime -= 1440
    }
    var hourAngle = trueSolarTime / 4.0 - 180.0;
    if (hourAngle < -180) {
        hourAngle += 360.0
    }
    var haRad = degToRad(hourAngle)
    var csz = Math.sin(degToRad(latitude)) * Math.sin(degToRad(theta)) + Math.cos(degToRad(latitude)) * Math.cos(degToRad(theta)) * Math.cos(haRad)
    if (csz > 1.0) {
        csz = 1.0
    } else if (csz < -1.0) {
        csz = -1.0
    }
    var zenith = radToDeg(Math.acos(csz))
    var azDenom = (Math.cos(degToRad(latitude)) * Math.sin(degToRad(zenith)))
    if (Math.abs(azDenom) > 0.001) {
        var azRad = ((Math.sin(degToRad(latitude)) * Math.cos(degToRad(zenith))) - Math.sin(degToRad(theta))) / azDenom
        if (Math.abs(azRad) > 1.0) {
            if (azRad < 0) {
                azRad = -1.0
            } else {
                azRad = 1.0
            }
        }
        var azimuth = 180.0 - radToDeg(Math.acos(azRad))
        if (hourAngle > 0.0) {
            azimuth = -azimuth
        }
    } else {
        if (latitude > 0.0) {
            var azimuth = 180.0
        } else {
            var azimuth = 0.0
        }
    }
    if (azimuth < 0.0) {
        azimuth += 360.0
    }
    var exoatmElevation = 90.0 - zenith

    // Atmospheric Refraction correction
    var refractionCorrection = calcRefraction(exoatmElevation)

    var solarZen = zenith - refractionCorrection;
    var elevation = 90.0 - solarZen

    return { "azimuth": azimuth, "elevation": elevation }
}

/**
 * Calculates the solar noon (local apparent noon) for a given Julian date, longitude, and timezone.
 *
 * @param {number} jd - The Julian date.
 * @param {number} longitude - The longitude in degrees (positive for east, negative for west).
 * @param {number} timezone - The timezone offset from UTC in hours.
 * @returns {number} The solar noon in local time (minutes from midnight).
 */
function calcSolNoon(jd, longitude, timezone) {
    var tnoon = calcTimeJulianCent(jd - longitude / 360.0)
    var eqTime = calcEquationOfTime(tnoon)
    var solNoonOffset = 720.0 - (longitude * 4) - eqTime // in minutes
    var newt = calcTimeJulianCent(jd - 0.5 + solNoonOffset / 1440.0)
    eqTime = calcEquationOfTime(newt)
    var solNoonLocal = 720 - (longitude * 4) - eqTime + (timezone * 60.0)// in minutes
    while (solNoonLocal < 0.0) {
        solNoonLocal += 1440.0;
    }
    while (solNoonLocal >= 1440.0) {
        solNoonLocal -= 1440.0;
    }

    return solNoonLocal
}

/**
 * Calculates the Universal Time Coordinated (UTC) time of sunrise or sunset for a given date and location.
 *
 * @param {boolean} rise - If 1, calculates the time of sunrise; if 0, calculates the time of sunset.
 * @param {number} JD - The Julian Date for which the calculation is to be made.
 * @param {number} latitude - The latitude of the location in degrees.
 * @param {number} longitude - The longitude of the location in degrees.
 * @param {number} horizon - The zenith angle of the horizon in degrees (e.g., 90.833 for the standard value of the sun's upper limb touching the horizon).
 * @returns {number} The time of sunrise or sunset in minutes from midnight UTC.
 */
function calcSunUTC(rise, JD, latitude, longitude, horizon) {
    var t = calcTimeJulianCent(JD);
    var eqTime = calcEquationOfTime(t);
    var solarDec = calcSunDec(t);
    var hourAngle = calcHourAngleSun(latitude, solarDec, horizon);
    if (!rise) hourAngle = -hourAngle;
    var delta = longitude + radToDeg(hourAngle);
    var timeUTC = 720 - (4.0 * delta) - eqTime;	// in minutes

    return timeUTC
}

// rise = 1 for sunrise, 0 for sunset
/**
 * Calculates the time and azimuth of sunrise or sunset.
 *
 * @param {boolean} rise - True for sunrise, false for sunset.
 * @param {number} JD - Julian Date.
 * @param {number} latitude - Latitude in degrees.
 * @param {number} longitude - Longitude in degrees.
 * @param {number} horizon - The zenith angle of the horizon in degrees
 * @param {number} timezone - Timezone offset from UTC in hours.
 * @returns {Object} An object containing:
 *   - {number} jday - Julian day of the event.
 *   - {number} timelocal - Local time of the event in minutes from midnight.
 *   - {number} azimuth - Azimuth of the sun at the event in degrees.
 */
function calcSun(rise, JD, latitude, longitude, horizon, timezone) {

    var timeUTC = calcSunUTC(rise, JD, latitude, longitude, horizon);
    var newTimeUTC = calcSunUTC(rise, JD + timeUTC / 1440.0, latitude, longitude, horizon);
    if (isNumber(newTimeUTC)) {
        var timeLocal = newTimeUTC + (timezone * 60.0)
        var riseT = calcTimeJulianCent(JD + newTimeUTC / 1440.0)
        var riseAzEl = calcAzEl(riseT, timeLocal, latitude, longitude, timezone)
        var azimuth = riseAzEl.azimuth
        var jday = JD
        if ((timeLocal < 0.0) || (timeLocal >= 1440.0)) {
            var increment = ((timeLocal < 0) ? 1 : -1)
            while ((timeLocal < 0.0) || (timeLocal >= 1440.0)) {
                timeLocal += increment * 1440.0
                jday -= increment
            }
        }

    } else { // no sunrise/set found

        var azimuth = -1.0
        var timeLocal = 0.0
        var doy = calcDoyFromJD(JD)
        if (((latitude > 66.4) && (doy > 79) && (doy < 267)) ||
            ((latitude < -66.4) && ((doy < 83) || (doy > 263)))) {
            //previous sunrise/next sunset
            jday = calcJDofNextPrev(!rise, rise, JD, latitude, longitude, timezone)
        } else {   //previous sunset/next sunrise
            jday = calcJDofNextPrev(rise, rise, JD, latitude, longitude, timezone)
        }
    }

    return { "jday": jday, "timelocal": timeLocal, "azimuth": azimuth }
}

/**
 * Calculates the Julian Date (JD) of the next or previous sunrise or sunset.
 *
 * @param {boolean} next - If true, calculates the next event; if false, calculates the previous event.
 * @param {boolean} rise - If true, calculates the sunrise; if false, calculates the sunset.
 * @param {number} JD - The initial Julian Date.
 * @param {number} latitude - The latitude of the observer in degrees.
 * @param {number} longitude - The longitude of the observer in degrees.
 * @param {number} horizon - The zenith angle of the horizon in degrees.
 * @param {number} tz - The time zone offset from UTC in hours.
 * @returns {number} - The Julian Date of the next or previous sunrise or sunset.
 */
function calcJDofNextPrev(next, rise, JD, latitude, longitude, horizon, tz) {

    var julianday = JD;
    var increment = ((next) ? 1.0 : -1.0);
    var time = calcSunUTC(rise, julianday, latitude, longitude, horizon);

    while (!isNumber(time)) {
        julianday += increment;
        time = calcSunUTC(rise, julianday, latitude, longitude, horizon);
    }
    var timeLocal = time + tz * 60.0
    while ((timeLocal < 0.0) || (timeLocal >= 1440.0)) {
        var incr = ((timeLocal < 0) ? 1 : -1)
        timeLocal += (incr * 1440.0)
        julianday -= incr
    }

    return julianday;
}

/******************************************************************/
/* Main functions that return a string in Stellarium date format. */
/******************************************************************/

/**
 * Calculates the local time of sunrise for a given Julian Date, latitude, longitude, and time zone.
 *
 * @param {number} JD - The Julian Date for which to calculate the sunrise time.
 * @param {number} lat - The observer latitude in degrees.
 * @param {number} lon - The observer longitude in degrees.
 * @param {number} tz - The time zone offset from UTC in hours.
 * @returns {string} The local time of sunrise (YYYY-MM-DDTHH:MM:SS).
 */
function sunRise(JD, lat, lon, tz) {
    t = calcSun(1, JD, lat, lon, 90.8333, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the local time of sunset for a given Julian Date, latitude, longitude, and time zone.
 *
 * @param {number} JD - The Julian Date for which to calculate the sunset time.
 * @param {number} lat - The latitude of the location in degrees.
 * @param {number} lon - The longitude of the location in degrees.
 * @param {number} tz - The time zone offset from UTC in hours.
 * @returns {string} The local time of sunset (YYYY-MM-DDTHH:MM:SS).
 */
function sunSet(JD, lat, lon, tz) {
    t = calcSun(0, JD, lat, lon, 90.8333, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the time of morning civil twilight.
 *
 * @param {number} JD - Julian Date.
 * @param {number} lat - Latitude in degrees.
 * @param {number} lon - Longitude in degrees.
 * @param {number} tz - Time zone offset from UTC in hours.
 * @returns {string} The local time of morning civil twilight (YYYY-MM-DDTHH:MM:SS).
 */
function civTwilightMorning(JD, lat, lon, tz) {
    t = calcSun(1, JD, lat, lon, 96, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the time of evening civil twilight.
 *
 * @param {number} JD - Julian Date.
 * @param {number} lat - Latitude in degrees.
 * @param {number} lon - Longitude in degrees.
 * @param {number} tz - Time zone offset from UTC in hours.
 * @returns {string} The local time of evening civil twilight (YYYY-MM-DDTHH:MM:SS).
 */
function civTwilightEvening(JD, lat, lon, tz) {
    t = calcSun(0, JD, lat, lon, 96, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the time of morning nautical twilight.
 *
 * @param {number} JD - Julian Date.
 * @param {number} lat - Latitude of the observer in degrees.
 * @param {number} lon - Longitude of the observer in degrees.
 * @param {number} tz - Time zone offset from UTC in hours.
 * @returns {string} The local time of morning nautical twilight (YYYY-MM-DDTHH:MM:SS).
 */
function nautTwilightMorning(JD, lat, lon, tz) {
    t = calcSun(1, JD, lat, lon, 102, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the time of evening nautical twilight.
 *
 * @param {number} JD - Julian Date.
 * @param {number} lat - Latitude in degrees.
 * @param {number} lon - Longitude in degrees.
 * @param {number} tz - Time zone offset from UTC in hours.
 * @returns {string} The local time of evening nautical twilight (YYYY-MM-DDTHH:MM:SS).
 */
function nautTwilightEvening(JD, lat, lon, tz) {
    t = calcSun(0, JD, lat, lon, 102, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the time of morning astronomical twilight.
 *
 * @param {number} JD - Julian Date.
 * @param {number} lat - Latitude in degrees.
 * @param {number} lon - Longitude in degrees.
 * @param {number} tz - Time zone offset from UTC in hours.
 * @returns {string} The local time of morning astronomical twilight (YYYY-MM-DDTHH:MM:SS).
 */
function astroTwilightMorning(JD, lat, lon, tz) {
    t = calcSun(1, JD, lat, lon, 108, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/**
 * Calculates the time of evening astronomical twilight.
 *
 * @param {number} JD - Julian Date.
 * @param {number} lat - Latitude in degrees.
 * @param {number} lon - Longitude in degrees.
 * @param {number} tz - Time zone offset from UTC in hours.
 * @returns {string} The local time of evening astronomical twilight (YYYY-MM-DDTHH:MM:SS).
 */
function astroTwilightEvening(JD, lat, lon, tz) {
    t = calcSun(0, JD, lat, lon, 108, tz);
    return timeDateString(t["jday"], t["timelocal"]);
}

/********************************************/
/* Some functions to deal with date strings */
/********************************************/

monthList = [
    { name: "January", numdays: 31, abbr: "Jan" },
    { name: "February", numdays: 28, abbr: "Feb" },
    { name: "March", numdays: 31, abbr: "Mar" },
    { name: "April", numdays: 30, abbr: "Apr" },
    { name: "May", numdays: 31, abbr: "May" },
    { name: "June", numdays: 30, abbr: "Jun" },
    { name: "July", numdays: 31, abbr: "Jul" },
    { name: "August", numdays: 31, abbr: "Aug" },
    { name: "September", numdays: 30, abbr: "Sep" },
    { name: "October", numdays: 31, abbr: "Oct" },
    { name: "November", numdays: 30, abbr: "Nov" },
    { name: "December", numdays: 31, abbr: "Dec" },
];

/**
 * Converts Julian Date and minutes into a formatted date-time string.
 *
 * @param {number} JD - The Julian Date.
 * @param {number} minutes - The number of minutes past midnight.
 * @returns {string} A formatted date-time string.
 */
function timeDateString(JD, minutes) {
    var date = calcDateFromJD(JD);
    var time = timeString(minutes, 2);
    return getDateString({
        year: date.year,
        month: date.month,
        day: Math.floor(date.day),
        hour: Math.floor(minutes / 60),
        minute: Math.floor(minutes % 60),
        second: 0
    });
}

/**
 * Converts a given number of minutes into a formatted time string.
 *
 * @param {number} minutes - The number of minutes to convert. Should be between 0 and 1440.
 * @param {number} flag - Determines the format of the output string.
 *                        If flag > 2, the output includes seconds.
 *                        If flag == 2 and seconds are 30 or more, the minute is rounded up.
 * @returns {string} The formatted time string in "HH:MM" or "HH:MM:SS" format, or "error" if the input is out of range.
 */
function timeString(minutes, flag) {
    if ((minutes >= 0) && (minutes < 1440)) {
        var floatHour = minutes / 60.0;
        var hour = Math.floor(floatHour);
        var floatMinute = 60.0 * (floatHour - Math.floor(floatHour));
        var minute = Math.floor(floatMinute);
        var floatSec = 60.0 * (floatMinute - Math.floor(floatMinute));
        var second = Math.floor(floatSec + 0.5);
        if (second > 59) {
            second = 0
            minute += 1
        }
        if ((flag == 2) && (second >= 30)) minute++;
        if (minute > 59) {
            minute = 0
            hour += 1
        }
        var output = zeroPad(hour, 2) + ":" + zeroPad(minute, 2);
        if (flag > 2) output = output + ":" + zeroPad(second, 2);
    } else {
        var output = "error"
    }

    return output;
}

/**
 * Pads a number with leading zeros until it reaches the specified number of digits.
 *
 * @param {number|string} n - The number to be padded.
 * @param {number} digits - The desired length of the resulting string.
 * @returns {string} The padded number as a string.
 */
function zeroPad(n, digits) {

    n = n.toString();
    while (n.length < digits) {
        n = '0' + n;
    }
    return n;
}

//--------------------------------------------------------------
/**
 * Converts a date object to an ISO 8601 formatted string.
 *
 * @param {Object} date - The date object to format.
 * @param {number} date.year - The year.
 * @param {number} date.month - The month (1-12).
 * @param {number} date.day - The day of the month (1-31).
 * @param {number} date.hour - The hour (0-23).
 * @param {number} date.minute - The minute (0-59).
 * @param {number} date.second - The second (0-59).
 * @returns {string} The formatted date string in ISO 8601 format.
 */
function getDateString(date) {

    var s = date.year
        + '-'
        + zeroPad(date.month, 2)
        + '-'
        + zeroPad(date.day, 2)
        + 'T'
        + zeroPad(date.hour, 2)
        + ':'
        + zeroPad(date.minute, 2)
        + ':'
        + zeroPad(date.second, 2)

    return s
}