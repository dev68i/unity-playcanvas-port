export const E       = globalThis.Math.E;
export const LN10    = globalThis.Math.LN10;
export const LN2     = globalThis.Math.LN2;
export const LOG2E   = globalThis.Math.LOG2E;
export const LOG10E  = globalThis.Math.LOG10E;
export const PI      = globalThis.Math.PI;
export const SQRT1_2 = globalThis.Math.SQRT1_2;
export const SQRT2   = globalThis.Math.SQRT2;
export const infinity = globalThis.Infinity;

export const positiveInfinity = Number.POSITIVE_INFINITY;
export const negativeInfinity = Number.NEGATIVE_INFINITY;

export const abs    = globalThis.Math.abs;
export const acos   = globalThis.Math.acos;
export const asin   = globalThis.Math.asin;
export const atan   = globalThis.Math.atan;
export const atan2  = globalThis.Math.atan2;
export const ceil   = globalThis.Math.ceil;
export const cos    = globalThis.Math.cos;
export const exp    = globalThis.Math.exp;
export const floor  = globalThis.Math.floor;
export const log    = globalThis.Math.log;
export const max    = globalThis.Math.max;
export const min    = globalThis.Math.min;
export const pow    = globalThis.Math.pow;
export const round  = globalThis.Math.round
export const sin    = globalThis.Math.sin;
export const sqrt   = globalThis.Math.sqrt;
export const tan    = globalThis.Math.tan;
export const random = globalThis.Math.random;

export const epsilon = Number.EPSILON;

export const kEpsilon = 0.00001;
export const kEpsilonNormalSqrt = 1e-15;

export const sqrt3Inv = 1.0 / sqrt(3.0);

/**
 * Degrees-to-radians conversion constant (Read Only).
 */
export const deg2Rad = PI / 180;

/**
 * Radians-to-degrees conversion constant (Read Only).
 */
export const rad2Deg = 360 / (PI * 2);

/**
 * Returns the sign of f.
 */
export const sign = (f: number) => (f >= 0) ? 1 : (-1);

/**
 * Loops the value t, so that it is never larger than length and never smaller than 0.
 * @param t 
 * @param length 
 * @returns 
 */
export const repeat = (t: number, length: number) => clamp(t - floor(t / length) * length, 0, length);

/**
 * Clamps the given value between the given minimum float and maximum float values.
 * Returns the given value if it is within the minimum and maximum range.
 * @param value The floating point value to restrict inside the range defined by the minimum and maximum values.
 * @param min The minimum floating point value to compare against.
 * @param max The maximum floating point value to compare against.
 * @returns The float result between the minimum and maximum values.
 */
export const clamp = (value: number, min: number, max: number) => value < min ? min : value > max ? max : value;

/**
 * Clamps value between 0 and 1 and returns value.
 * @returns 
 */
export const clamp01 = (value: number) => value < 0 ? 0 : value > 1 ? 1 : value;

/**
 * Linearly interpolates between a and b by t.
 * @param a The start value.
 * @param b The end value.
 * @param t The interpolation value between the two floats.
 * @returns The interpolated float result between the two float values.
 */
export const lerp = (a: number, b: number, t: number) => a + (b - a) * clamp01(t);

/**
 * Linearly interpolates between a and b by t with no limit to t.
 * @param a The start value.
 * @param b The end value.
 * @param t The interpolation between the two floats.
 * @returns The float value as a result from the linear interpolation.
 */
export const lerpUnclamped = (a: number, b: number, t: number) => a + (b - a) * t;

/**
 * Same as Lerp but makes sure the values interpolate correctly when they wrap around 360 degrees.
 * @param a The start angle. A float expressed in degrees.
 * @param b The end angle. A float expressed in degrees.
 * @param t The interpolation value between the start and end angles. This value is clamped to the range [0, 1].
 * @returns Returns the interpolated float result between angle a and angle b, based on the interpolation value t.
 */
export const lerpAngle = (a: number, b: number, t: number) => {

    let num = repeat(b - a, 360);

    if (num > 180)
    {
        num -= 360;
    }

    return a + num * clamp01(t);
}

/**
 * Moves a value current towards target.
 * @param current The current value.
 * @param target The value to move towards.
 * @param maxDelta The maximum change that should be applied to the value.
 * @returns 
 */
export const moveTowards = (current: number, target: number, maxDelta: number) => {

    if (abs(target - current) <= maxDelta){
        return target;
    }

    return current + sign(target - current) * maxDelta;
}

/**
 * Same as MoveTowards but makes sure the values interpolate correctly when they wrap around 360 degrees.
 */
export const moveTowardsAngle = (current: number, target: number, maxDelta: number) => {

    const num = deltaAngle(current, target);

    if (0 - maxDelta < num && num < maxDelta) {
        return target;
    }

    target = current + num;
    return moveTowards(current, target, maxDelta);
}

/**
 * Interpolates between min and max with smoothing at the limits.
 */
export const smoothStep = (from: number, to: number, t: number) => {
    t = clamp01(t);
    t = -2 * t * t * t + 3 * t * t;
    return to * t + from * (1 - t);
}

export const gamma = (value: number, absmax: number, gamma: number) => {

    const flag = value < 0;
    const num = abs(value);

    if (num > absmax) {
        return flag ? (0 - num) : num;
    }

    const num2 = pow(num / absmax, gamma) * absmax;
    return flag ? (0 - num2) : num2;
}

/**
 * Compares two floating point values and returns true if they are similar.
 */
export const approximately = (a: number, b: number) => {
    return abs(b - a) < max(1E-06 * max(abs(a), abs(b)), epsilon * 8);
}

/**
 * Gradually changes a value towards a desired goal over time.
 */
export const smoothDamp = (
    current: number,
    target: number,
    refCurrentVelocity: { currentVelocity: number},
    smoothTime: number,
    maxSpeed: number,
    deltaTime: number
) => {

    // Based on Game Programming Gems 4 Chapter 1.10
    smoothTime = max(0.0001, smoothTime);
    const omega = 2 / smoothTime;

    const x = omega * deltaTime;
    const exp = 1 / (1 + x + 0.48 * x * x + 0.235 * x * x * x);
    let change = current - target;
    const originalTo = target;

    // Clamp maximum speed
    const maxChange = maxSpeed * smoothTime;
    change = clamp(change, -maxChange, maxChange);
    target = current - change;

    const temp = (refCurrentVelocity.currentVelocity + omega * change) * deltaTime;
    refCurrentVelocity.currentVelocity = (refCurrentVelocity.currentVelocity - omega * temp) * exp;

    let output = target + (change + temp) * exp;

    // Prevent overshooting
    if (originalTo - current > 0.0 === output > originalTo) {
        output = originalTo;
        refCurrentVelocity.currentVelocity = (output - originalTo) / deltaTime;
    }

    return output;
}

/**
 * PingPong returns a value that will increment and decrement between the value 0 and length.
 */
export const pingPong = (t: number, length: number) => {
    t = repeat(t, length * 2);
    return length - abs(t - length);
}

/**
 * Determines where a value lies between two points.
 */
export const inverseLerp = (a: number, b: number, value: number) => a != b ? clamp01((value - a) / (b - a)) : 0;

/**
 * Calculates the shortest difference between two given angles given in degrees.
 */
export const deltaAngle = (current: number, target: number) => {
    let num = repeat(target - current, 360);
    if (num > 180)
    {
        num -= 360;
    }
    return num;
}

/**
 * Rounds the value to the specified value in an area with a given radius from 0
 */
export const roundAtArea = (current: number, raduis: number, target: number) => {

    if (current > raduis ||
        current < -raduis) {
        return current;
    }

    return target;
}

/**
 * Is floatA equal to zero? Takes floating point inaccuracy into account, by using E^-3.
 * @param floatA 
 * @returns 
 */
export const isEqualToZero = (floatA: number, e: number = 1e-3): boolean => {
    return abs(floatA) < e;
}
		
/**
 * Is floatA not equal to zero? Takes floating point inaccuracy into account, by using E^-3.
 */
export const notEqualToZero = (floatA: number, e: number = 1e-3): boolean => {
    return abs(floatA) > e;
}

export const Math = {

    E,
    LN10,
    LN2,
    LOG2E,
    LOG10E,
    PI,
    SQRT1_2,
    SQRT2,

    abs,
    acos,
    asin,
    atan,
    atan2,
    ceil,
    cos,
    exp,
    floor,
    log,
    max,
    min,
    pow,
    round,
    sin,
    sqrt,
    tan,
    random,

    infinity,
    positiveInfinity,
    negativeInfinity,
    sqrt3Inv,
    epsilon,
    deg2Rad,
    rad2Deg,
    sign,
    repeat,
    clamp,
    clamp01,
    lerp,
    lerpUnclamped,
    lerpAngle,
    moveTowards,
    moveTowardsAngle,
    smoothStep,
    gamma,
    approximately,
    smoothDamp,
    pingPong,
    inverseLerp,
    deltaAngle,
    roundAtArea,

    isEqualToZero,
    notEqualToZero,
}

export default Math;
