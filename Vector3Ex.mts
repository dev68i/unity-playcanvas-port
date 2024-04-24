import { type IVector3, Quaternion, Vector3 } from "./Types.mjs";
import Math, { kEpsilonNormalSqrt } from './Math.mjs';
import { UnityVector3 } from "./Vector3Dfl.mjs";
import QuaternionExBase from "./QuaternionExBase.mjs";

export type IMutableVector3 = IVector3;

/**
 * Sets the specified 3-dimensional vector to the supplied numerical values.
 */
export function set<T extends IMutableVector3>(to: T, x: number, y: number, z: number): T {
    to.x = x;
    to.y = y;
    to.z = z;
    return to;
}

/**
 * Copies the contents of a source 3-dimensional vector to a destination 3-dimensional vector.
 */
export function assign<T extends IMutableVector3>(to: T, from: Readonly<IVector3>): T {
    to.x = from.x;
    to.y = from.y;
    to.z = from.z;
    return to;
}

/**
 * Returns an identical copy of the specified 3-dimensional vector.
 */
export function clone(from: Readonly<IVector3>): Vector3 {
    return new Vector3(from.x, from.y, from.z);
}

/**
 * Returns a vector with the same direction as the given vector, but with a length of 1.
 */
export function normalize(value: Readonly<IVector3>): Vector3;
export function normalize<TOut extends IMutableVector3>(value: Readonly<IVector3>, out: TOut): TOut;
export function normalize(value: Readonly<IVector3>, out = new Vector3()): Vector3 {

    const ls = value.x * value.x + value.y * value.y + value.z * value.z;
    const length = Math.sqrt(ls);

    if(length !== 0) {
        return set(
            out,
            value.x / length,
            value.y / length,
            value.z / length
        );
    }

    return assign(out, value);
}

/**
 * Update the vector with the same direction as the given vector, but with a length of 1.
 */
export function normalizeRef<T extends IMutableVector3>(refValue: T): T {

    const ls = refValue.x * refValue.x + refValue.y * refValue.y + refValue.z * refValue.z;
    const length = Math.sqrt(ls);

    if(length !== 0) {
        refValue.x /= length;
        refValue.y /= length;
        refValue.z /= length;
        return refValue;
    }
    
    return refValue;
}

/**
 * Returns a vector with the same direction as the given vector, but with a set magnitude.
 * If value eq zero return zero vector
 */
export function setMagnitude(value: Readonly<IVector3>, length: number): Vector3;
export function setMagnitude<TOut extends IMutableVector3>(value: Readonly<IVector3>, length: number, out: TOut): TOut;
export function setMagnitude(value: Readonly<IVector3>, length: number, out: IVector3 = new Vector3()) {

    if ((value.x === 0 && value.z === 0 && value.z === 0) || length === 0) {
        return set(out, 0, 0, 0);
    }

    const ls = value.x * value.x + value.y * value.y + value.z * value.z;
    const vecLength = Math.sqrt(ls);

    if(length !== 0) {
        return set(
            out,
            value.x / vecLength * length,
            value.y / vecLength * length,
            value.z / vecLength * length
        );
    }

    return assign(out, value);
}

/**
 * Returns a vector with the same direction as the given vector, but with a set magnitude.
 */
export function setMagnitudeRef<T extends IMutableVector3>(refValue: T, length: number): T {
    return setMagnitude(refValue, length, refValue);
}

/**
 * Computes the cross product of two vectors.
 */
export function cross(vector1: Readonly<IVector3>, vector2: Readonly<IVector3>): Vector3;
export function cross<TOut extends IMutableVector3>(vector1: Readonly<IVector3>, vector2: Readonly<IVector3>, out: TOut): TOut;
export function cross(vector1: Readonly<IVector3>, vector2: Readonly<IVector3>, out = new Vector3()) {
    return set(
        out,
        vector1.y * vector2.z - vector1.z * vector2.y,
        vector1.z * vector2.x - vector1.x * vector2.z,
        vector1.x * vector2.y - vector1.y * vector2.x
    );
}

/**
 * Computes the cross product of two vectors.
 */
export function crossRef<T extends IMutableVector3>(refVector1: T, vector2: Readonly<IVector3>): T {
    return cross(refVector1, vector2, refVector1);
}

/**
 * Computes the cross product of two vectors.
 */
export function crossRef2<T extends IMutableVector3>(vector1: Readonly<IVector3>, refVector2: T): T {
    return cross(vector1, refVector2, refVector2);
}

/**
 * Returns the reflection of a vector off a surface that has the specified normal.
 */
export function reflect(vector: Readonly<IVector3>, normal: Readonly<IVector3>): Vector3;
export function reflect<TOut extends IMutableVector3>(vector: Readonly<IVector3>, normal: Readonly<IVector3>, out: TOut): TOut;
export function reflect(vector: Readonly<IVector3>, normal: Readonly<IVector3>, out = new Vector3()) {
    const dot   = vector.x * normal.x + vector.y * normal.y + vector.z * normal.z;
    const tempX = normal.x * dot * 2;
    const tempY = normal.y * dot * 2;
    const tempZ = normal.z * dot * 2;
    return set(
        out,
        vector.x - tempX,
        vector.y - tempY,
        vector.z - tempZ
    );
}

/**
 * Returns the reflection of a vector off a surface that has the specified normal.
 */
export function reflectRef<T extends IMutableVector3>(refVector: T, normal: Readonly<IVector3>): T {
    return reflect(refVector, normal, refVector);
}

/**
 * Returns the reflection of a vector off a surface that has the specified normal.
 */
export function reflectRef2<T extends IMutableVector3>(vector: Readonly<IVector3>, refNormal: T): T {
    return reflect(vector, refNormal, refNormal);
}

/**
 * Returns the Euclidean distance between the two given points.
 */
export function distance(value1: Readonly<IVector3>, value2: Readonly<IVector3>) {
    const dx = value1.x - value2.x;
    const dy = value1.y - value2.y;
    const dz = value1.z - value2.z;
    const ls = dx * dx + dy * dy + dz * dz;
    return Math.sqrt(ls);
}

/**
 * Returns a copy of vector with its magnitude clamped to maxLength.
 */
export function clampMagnitude(vector: Readonly<IVector3>, maxLength: number): Vector3;
export function clampMagnitude<TOut extends IMutableVector3>(vector: Readonly<IVector3>, maxLength: number, out: TOut): TOut;
export function clampMagnitude(vector: Readonly<IVector3>, maxLength: number, out = new Vector3()) {
    
    const num = sqrMagnitude(vector);

    if (num > maxLength * maxLength) {

        const num2 = Math.sqrt(num);
        const num3 = vector.x / num2;
        const num4 = vector.y / num2;
        const num5 = vector.z / num2;

        return set(out, num3 * maxLength, num4 * maxLength, num5 * maxLength);
    }

    return assign(out, vector);
}

/**
 * Returns a copy of vector with its magnitude clamped to maxLength.
 */
export function clampMagnitudeRef<T extends IMutableVector3>(refVector: T, maxLength: number): T {
    return clampMagnitude(refVector, maxLength, refVector);
}

/**
 * Returns the Euclidean distance squared between the two given points.
 */
export function distanceSquared(value1: Readonly<Vector3>, value2: Readonly<Vector3>) {
    const dx = value1.x - value2.x;
    const dy = value1.y - value2.y;
    const dz = value1.z - value2.z;
    return dx * dx + dy * dy + dz * dz;
}

/**
 * Returns the dot product of two vectors.
 */
export function dot(vector1: Readonly<IVector3>, vector2: Readonly<IVector3>) {
    return  vector1.x * vector2.x +
            vector1.y * vector2.y +
            vector1.z * vector2.z;
}

/**
 * Returns a vector whose elements are the minimum of each of the pairs of elements in the two source vectors.
 */
export function min(value1: Readonly<Vector3>, value2: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        (value1.x < value2.x) ? value1.x : value2.x,
        (value1.y < value2.y) ? value1.y : value2.y,
        (value1.z < value2.z) ? value1.z : value2.z
    );
}

/**
 * Returns a vector whose elements are the minimum of each of the pairs of elements in the two source vectors.
 */
export function minRef(refValue1: Vector3, value2: Readonly<Vector3>) {
    return min(refValue1, value2, refValue1);
}

/**
 * Returns a vector whose elements are the maximum of each of the pairs of elements in the two source vectors.
 */
export function max(value1: Readonly<Vector3>, value2: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        (value1.x > value2.x) ? value1.x : value2.x,
        (value1.y > value2.y) ? value1.y : value2.y,
        (value1.z > value2.z) ? value1.z : value2.z
    );
}

/**
 * Returns a vector whose elements are the maximum of each of the pairs of elements in the two source vectors.
 */
export function maxRef(refValue1: Vector3, value2: Readonly<Vector3>) {
    return max(refValue1, value2, refValue1);
}

/**
 * Returns a vector whose elements are the maximum of each of the pairs of elements in the two source vectors.
 */
export function maxRef2(value1: Readonly<Vector3>, refValue2: Vector3) {
    return max(value1, refValue2, refValue2);
}

/**
 * Returns a vector whose elements are the absolute values of each of the source vector's elements.
 */
export function abs(value: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        Math.abs(value.x),
        Math.abs(value.y),
        Math.abs(value.z)
    );
}

/**
 * Returns a vector whose elements are the absolute values of each of the source vector's elements.
 */
export function absRef(refValue: Vector3) {
    return abs(refValue, refValue);
}

/**
 * Restricts a vector between a min and max value.
 */
export function clamp(value1: Readonly<Vector3>, min: Readonly<Vector3>, max: Readonly<Vector3>, out = new Vector3()) {

    // This compare order is very important!!!
    // We must follow HLSL behavior in the case user specified min value is bigger than max value.

    let x = value1.x;
    x = (x > max.x) ? max.x : x;
    x = (x < min.x) ? min.x : x;

    let y = value1.y;
    y = (y > max.y) ? max.y : y;
    y = (y < min.y) ? min.y : y;

    let z = value1.z;
    z = (z > max.z) ? max.z : z;
    z = (z < min.z) ? min.z : z;

    return set(out, x, y, z);
}

/**
 * Restricts a vector between a min and max value.
 */
export function clampRef(refValue1: Vector3, min: Readonly<Vector3>, max: Readonly<Vector3>) {
    return clamp(refValue1, min, max, refValue1);
}

/**
 * Restricts a vector between a min and max value.
 */
export function clampRef2(value: Readonly<Vector3>, refMin: Vector3, max: Readonly<Vector3>) {
    return clamp(value, refMin, max, refMin);
}

/**
 * Restricts a vector between a min and max value.
 */
export function clampRef3(value: Readonly<Vector3>, min: Readonly<Vector3>, refMax: Vector3) {
    return clamp(value, min, refMax, refMax);
}

/**
 * Linearly interpolates between two vectors based on the given weighting.
 */
export function lerpUnclamped(value1: Readonly<Vector3>, value2: Readonly<Vector3>, amount: number, out = new Vector3()) {
    return set(
        out,
        value1.x + (value2.x - value1.x) * amount,
        value1.y + (value2.y - value1.y) * amount,
        value1.z + (value2.z - value1.z) * amount
    );
}

/**
 * Linearly interpolates between two vectors based on the given weighting.
 */
export function lerpUnclampedRef(refValue1: Vector3, value2: Readonly<Vector3>, amount: number) {
    return lerpUnclamped(refValue1, value2, amount, refValue1);
}

/**
 * Linearly interpolates between two vectors based on the given weighting.
 */
export function lerpUnclampedRef2(value1: Readonly<Vector3>, refValue2: Vector3, amount: number) {
    return lerpUnclamped(value1, refValue2, amount, refValue2);
}

/**
 * Linearly interpolates between two vectors based on the given weighting [0..1].
 */
export function lerp(start: Readonly<Vector3>, end: Readonly<Vector3>, percent: number, out = new Vector3()) {
    percent = Math.clamp01(percent);
    return lerpUnclamped(start, end, percent, out);
}

/**
 * Linearly interpolates between two vectors based on the given weighting [0..1].
 */
export function lerpRef(refStart: Vector3, end: Readonly<Vector3>, percent: number) {
    return lerp(refStart, end, percent, refStart);
}

/**
 * Linearly interpolates between two vectors based on the given weighting [0..1].
 */
export function lerpRef2(start: Readonly<Vector3>, refEnd: Vector3, percent: number) {
    return lerp(start, refEnd, percent, refEnd);
}

export function findOrthogonal(v: Readonly<Vector3>, out = new Vector3()) {

    if (v.x >= Math.sqrt3Inv) {
        return set(out, v.y, -v.x, 0.0);
    }

    return set(out, 0.0, v.z, -v.y);
}

export function findOrthogonalRef(refV: Vector3) {
    return findOrthogonal(refV, refV);
}

/**
 * Multiplies a Quaternion by a Vector3 value.
 */
export function multiplyQV3(rotation: Readonly<Quaternion>, point: Readonly<Vector3>, out = new Vector3()) {

    const x = point.x;
    const y = point.y;
    const z = point.z;
    const qx = rotation.x;
    const qy = rotation.y;
    const qz = rotation.z;
    const qw = rotation.w;

    // q * v
    const ix =  qw * x + qy * z - qz * y;
    const iy =  qw * y + qz * x - qx * z;
    const iz =  qw * z + qx * y - qy * x;
    const iw = -qx * x - qy * y - qz * z;

    return set(
        out,
        ix * qw + iw * -qx + iy * -qz - iz * -qy,
        iy * qw + iw * -qy + iz * -qx - ix * -qz,
        iz * qw + iw * -qz + ix * -qy - iy * -qx,
    );

    /*
    const x  = rotation.x * 2;
    const y  = rotation.y * 2;
    const z  = rotation.z * 2;
    const xx = rotation.x * x;
    const yy = rotation.y * y;
    const zz = rotation.z * z;
    const xy = rotation.x * y;
    const xz = rotation.x * z;
    const yz = rotation.y * z;
    const wx = rotation.w * x;
    const wy = rotation.w * y;
    const wz = rotation.w * z;

    return set(
        out,
        (1 - (yy + zz)) * point.x + (xy - wz) * point.y + (xz + wy) * point.z,
        (xy + wz) * point.x + (1 - (xx + zz)) * point.y + (yz - wx) * point.z,
        (xz - wy) * point.x + (yz + wx) * point.y + (1 - (xx + yy)) * point.z
    );
    */
}

/**
 * Multiplies a Quaternion by a Vector3 value.
 */
export function multiplyQV3Ref2(rotation: Readonly<Quaternion>, refPoint: Vector3) {
    return multiplyQV3(rotation, refPoint, refPoint);
}

export function fastOrthogonal(v: Readonly<Vector3>, normalize = true, out = new Vector3()) {

    let sqr = v.x * v.x + v.y * v.y;

    // (0,0,1) x (x,y,z)
    if(sqr > 0) {
      const im = normalize ? 1 / Math.sqrt(sqr) : 1;
      return set(out, -v.y * im, v.x * im, 0);
    }
    
    // (1,0,0) x (x,y,z)
    sqr = v.y * v.y + v.z * v.z;
    const im = normalize? 1 / Math.sqrt(sqr) : 1;
    return set(out, 0, -v.z * im, v.y * im);
}

export function fastOrthogonalRef(refV: Vector3, normalize = true) {
    return fastOrthogonal(refV, normalize, refV);
}

const tmpSlerpUnclampedV3 = new Vector3();

/**
 * Spherically interpolates between two vectors.
 */
export function slerpUnclamped(from: Vector3, to: Vector3, t: number, out = new Vector3()) {

    const len1 = magnitude(from);
	const len2 = magnitude(to);

    const v1 = assign(tmpSlerpUnclampedV3, from);
	const v2 = assign(out, to);

    divideScalarRef(v2, len2);
    divideScalarRef(v1, len1);

	const len 	= (len2 - len1) * t + len1;
	const cosom = dot(v1, v2);
	
    let scale0 = 1 - t;
    let scale1 = t;

	if (1 - cosom > 1e-6) {

		const omega = Math.acos(cosom);
		const sinom = Math.sin(omega);

		scale0 	= Math.sin((1 - t) * omega) / sinom;
		scale1 	= Math.sin(t * omega) / sinom;
    }

    multiplyScalarRef(v1, scale0);
    multiplyScalarRef(v2, scale1);
    addRef(v2, v1);
    multiplyScalarRef(v2, len);

	return out;
}

/**
 * Spherically interpolates between two vectors.
 */
export function slerpUnclampedRef(refFrom: Vector3, to: Vector3, t: number) {
    return slerpUnclamped(refFrom, to, t, refFrom);
}

/**
 * Spherically interpolates between two vectors.
 */
export function slerpUnclampedRef2(from: Vector3, refTo: Vector3, t: number) {
    return slerpUnclamped(from, refTo, t, refTo);
}

/**
 * Spherically interpolates between two vectors.
 */
export function slerp(a: Vector3, b: Vector3, t: number, out = new Vector3()) {
	if (t <= 0) return assign(out, a);
    if (t >= 1) return assign(out, b);
    return slerpUnclamped(a, b, t, out);
}

/**
 * Spherically interpolates between two vectors.
 */
export function slerpRef(refA: Vector3, b: Vector3, t: number) {
    return slerp(refA, b, t, refA);
}

/**
 * Spherically interpolates between two vectors.
 */
export function slerpRef2(a: Vector3, refB: Vector3, t: number) {
    return slerp(a, refB, t, refB);
}

const tmpSlerp2UnclampedV3 = new Vector3();

/**
 * Spherically interpolates between two vectors. (Version 2)
 */
export function slerp2Unclamped(from: Vector3, to: Vector3, t: number, out = new Vector3()) {

    const len1 = magnitude(from);
    const len2 = magnitude(to);
    
    const v1 = assign(tmpSlerp2UnclampedV3, from);
    const v2 = assign(out, to);

    divideScalarRef(v2, len2);
    divideScalarRef(v1, len1);

    const omega = dot(v1, v2);
    const len 	= (len2 - len1) * t + len1;
    const theta = Math.acos(omega) * t;
    const tmp = multiplyScalar(v1, omega);

    subtractRef(v2, tmp);
    normalizeRef(v2);
    multiplyScalarRef(v2, Math.sin(theta));
    multiplyScalarRef(v1, Math.cos(theta));
    addRef(v2, v1);
    normalizeRef(v2);
    multiplyScalarRef(v2, len);

    return out;
}

/**
 * Spherically interpolates between two vectors. (Version 2)
 */
export function slerp2(a: Vector3, b: Vector3, t: number, out = new Vector3()) {
	if (t <= 0) return assign(out, a);
    if (t >= 1) return assign(out, b);
    return slerp2Unclamped(a, b, t, out);
}

/**
 * Returns a boolean indicating whether the two given vectors are equal.
 */
export function equals(left: Vector3, right: Vector3) { 
    return  left.x === right.x &&
            left.y === right.y &&
            left.z === right.z;
}

/**
 * Returns a boolean indicating whether the two given vectors are not equal.
 */
export function notEquals(left: Vector3, right: Vector3) {
    return  left.x !== right.x ||
            left.y !== right.y ||
            left.z !== right.z;
}

/**
 * Returns a vector whose elements are the square root of each of the source vector's elements.
 */
export function squareRoot(value: Vector3, out = new Vector3()) {
    return set(
        out,
        Math.sqrt(value.x),
        Math.sqrt(value.y),
        Math.sqrt(value.z)
    );
}

/**
 * Returns a vector whose elements are the square root of each of the source vector's elements.
 */
export function squareRootRef(refValue: Vector3) {
    return squareRoot(refValue, refValue);
}

/**
 * Adds two vectors.
 */
export function add(left: Vector3, right: Vector3, out = new Vector3()) {
    return set(
        out,
        left.x + right.x,
        left.y + right.y,
        left.z + right.z
    );
}

/**
 * Adds two vectors.
 */
export function addRef(refLeft: Vector3, right: Readonly<Vector3>) {
    return add(refLeft, right, refLeft);
}

/**
 * Adds two vectors.
 */
export function addRef2(left: Readonly<Vector3>, refRight: Vector3) {
    return add(left, refRight, refRight);
}

/**
 * Subtracts the second vector from the first.
 */
export function subtract(left: Readonly<Vector3>, right: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        left.x - right.x,
        left.y - right.y,
        left.z - right.z
    );
}

/**
 * Subtracts the second vector from the first.
 */
export function subtractRef(refLeft: Vector3, right: Readonly<Vector3>) {
    return subtract(refLeft, right, refLeft);
}

/**
 * Subtracts the second vector from the first.
 */
export function subtractRef2(left: Readonly<Vector3>, refRight: Vector3) {
    return subtract(left, refRight, refRight);
}

/**
 * Multiplies two vectors together.
 */
export function multiply(left: Readonly<Vector3>, right: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        left.x * right.x,
        left.y * right.y,
        left.z * right.z
    );
}

/**
 * Multiplies two vectors together.
 */
export function multiplyRef(refLeft: Vector3, right: Readonly<Vector3>) {
    return multiply(refLeft, right, refLeft);
}

/**
 * Multiplies two vectors together.
 */
export function multiplyRef2(left: Readonly<Vector3>, refRight: Vector3) {
    return multiply(left, refRight, refRight);
}

/**
 * Multiplies a vector by the given scalar.
 */
export function multiplyScalar(left: Readonly<Vector3>, right: number, out = new Vector3()) {
    return set(
        out,
        left.x * right,
        left.y * right,
        left.z * right
    );
}

/**
 * Multiplies a vector by the given scalar.
 */
export function multiplyScalarRef(refLeft: Vector3, right: number) {
    return multiplyScalar(refLeft, right, refLeft);
}

/**
 * Multiplies a vector by the given scalar.
 */
export function multiplyScalarL(left: number, right: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        left * right.x,
        left * right.y,
        left * right.z
    );
}

/**
 * Multiplies a vector by the given scalar.
 */
export function multiplyScalarLRef2(left: number, refRight: Vector3) {
    return multiplyScalarL(left, refRight, refRight);
}

/**
 * Divides the first vector by the second.
 */
export function divide(left: Vector3, right: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        left.x / right.x,
        left.y / right.y,
        left.z / right.z
    );
}

/**
 * Divides the first vector by the second.
 */
export function divideRef(refLeft: Vector3, right: Readonly<Vector3>) {
    return divide(refLeft, right, refLeft);
}

/**
 * Divides the first vector by the second.
 */
export function divideRef2(left: Readonly<Vector3>, refRight: Vector3) {
    return divide(left, refRight, refRight);
}

/**
 * Divides the vector by the given scalar.
 */
export function divideScalar(left: Readonly<Vector3>, right: number, out = new Vector3()) {
    return set(
        out,
        left.x / right,
        left.y / right,
        left.z / right
    );
}

/**
 * Divides the vector by the given scalar.
 */
export function divideScalarRef(refLeft: Vector3, right: number) {
    return divideScalar(refLeft, right, refLeft);
}

/**
 * Divides the scalar by the given vector.
 */
export function divideScalarL(left: number, right: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        left / right.x,
        left / right.y,
        left / right.z
    );
}

/**
 * Divides the scalar by the given vector.
 */
export function divideScalarLRef2(left: number, refRight: Vector3) {
    return divideScalarL(left, refRight, refRight);
}

/**
 * Negates a given vector.
 */
export function negate(value: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        -value.x,
        -value.y,
        -value.z
    );
}

/**
 * Negates a given vector.
 */
export function negateRef(refValue: Vector3) {
    return negate(refValue, refValue);
}

/**
 * Makes vectors normalized and orthogonal to each other.<br></br>
 * <b>Vectors modify by ref!!!</b>
 */
export function orthoNormalize(refNormal: Vector3, refTangent: Vector3) {

    normalizeRef(refNormal);

    const dotValue = dot(refTangent, refNormal);
    const vecValue = multiplyScalar(refNormal, dotValue);

    subtractRef(refTangent, vecValue);
    normalizeRef(refTangent);

    return {
        refNormal,
        refTangent
    }
}

export function angleRad(from: Readonly<IVector3>, to: Readonly<IVector3>) {

    // sqrt(a) * sqrt(b) = sqrt(a * b) -- valid for real numbers
    const denominator = Math.sqrt(sqrMagnitude(from) * sqrMagnitude(to));

    if (denominator < kEpsilonNormalSqrt) {
        return 0;
    }

    const fromToDot = dot(from, to);
    const dotValue = Math.clamp(fromToDot / denominator, -1, 1);

    return Math.acos(dotValue);
}

export function angle(from: Readonly<IVector3>, to: Readonly<IVector3>) {
    return angleRad(from, to) * Math.rad2Deg;
}

/**
 * Gets the closest direction axis to a vector. Input vector must be normalized!
 */
export function axis(v: Readonly<Vector3>, out = new Vector3()) {

    let neg = false;

    const x = dot(v, UnityVector3.right);
    const y = dot(v, UnityVector3.up);

    if (x < 0) neg = true;

    let absDot = Math.abs(y);
    let maxAbsDot = Math.abs(x);
    let closest = UnityVector3.right;

    if (absDot > maxAbsDot) {
        maxAbsDot = absDot;
        closest = UnityVector3.up;
        neg = y < 0;
    }

    const z = dot(v, UnityVector3.forward);
    absDot = Math.abs(z);

    if (absDot > maxAbsDot) {
        closest = UnityVector3.forward;
        neg = z < 0;
    }

    if (neg) {
        return set(out, -closest.x, -closest.y, -closest.z);
    }

    return assign(out, closest);
}

/**
 * Multiplies two vectors component-wise.
 */
export function scale(a: Readonly<Vector3>, b: Readonly<Vector3>, out = new Vector3()) {
    return set(
        out,
        a.x * b.x,
        a.y * b.y,
        a.z * b.z
    );
}

/**
 * Multiplies two vectors component-wise.
 */
export function scaleRef(refA: Vector3, b: Readonly<Vector3>) {
    return scale(refA, b, refA);
}

/**
 * Multiplies two vectors component-wise.
 */
export function scaleRef2(a: Readonly<Vector3>, refB: Vector3) {
    return scale(a, refB, refB);
}

/**
 * Projects a vector onto another vector.
 */
export function project(vector: Readonly<Vector3>, onNormal: Readonly<Vector3>, out = new Vector3()) {

    const num = dot(onNormal, onNormal);

    if (num < Math.epsilon) {
        return set(out, 0, 0, 0);
    }

    const num2  = dot(vector, onNormal);
    const right = num2 / num;

    return set(
        out,
        onNormal.x * right,
        onNormal.y * right,
        onNormal.z * right
    );
}

/**
 * Projects a vector onto another vector.
 */
export function projectRef(refVector: Vector3, onNormal: Readonly<Vector3>) {
    return project(refVector, onNormal, refVector);
}

/**
 * Projects a vector onto another vector.
 */
export function projectRef2(vector: Readonly<Vector3>, refOnNormal: Vector3) {
    return project(vector, refOnNormal, refOnNormal);
}

/**
 * Projects a vector onto a plane defined by a normal orthogonal to the plane.
 * @param vector The location of the vector above the plane.
 * @param planeNormal The direction from the vector towards the plane.
 */
export function projectOnPlane(vector: Readonly<Vector3>, planeNormal: Readonly<Vector3>, out: Vector3 = new Vector3()) {

    const num = dot(planeNormal, planeNormal);

    if (num < Math.epsilon) {
        return assign(out, vector);
    }

    const num2 = dot(vector, planeNormal);
    const right = num2 / num;
    
    return set(
        out,
        vector.x - planeNormal.x * right,
        vector.y - planeNormal.y * right,
        vector.z - planeNormal.z * right
    );
}

/**
 * Projects a vector onto a plane defined by a normal orthogonal to the plane.
 * @param refVector The location of the vector above the plane.
 * @param planeNormal The direction from the vector towards the plane.
 */
export function projectOnPlaneRef(refVector: Readonly<Vector3>, planeNormal: Readonly<Vector3>) {
    return projectOnPlane(refVector, planeNormal, refVector);
}

/**
 * Projects a vector onto a plane defined by a normal orthogonal to the plane.
 * @param vector The location of the vector above the plane.
 * @param refPlaneNormal The direction from the vector towards the plane.
 */
export function projectOnPlaneRef2(vector: Readonly<Vector3>, refPlaneNormal: Vector3) {
    return projectOnPlane(vector, refPlaneNormal, refPlaneNormal);
}

export function magnitude(vector: Readonly<IVector3>) {
    return Math.sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

export function sqrMagnitude(vector: Readonly<IVector3>) {
    return vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
}

/**
 * Calculates the signed angle between vectors from and to in relation to axis.
 * @param from The vector from which the angular difference is measured.
 * @param to The vector to which the angular difference is measured.
 * @param axis A vector around which the other vectors are rotated.
 * @returns Returns the signed angle between from and to in degrees.
*/
export function signedAngle(from: Readonly<IVector3>, to: Readonly<IVector3>, axis: Readonly<Vector3>) {
    const num  = angle(from, to);
    const num2 = from.y * to.z - from.z * to.y;
    const num3 = from.z * to.x - from.x * to.z;
    const num4 = from.x * to.y - from.y * to.x;
    const num5 = Math.sign(axis.x * num2 + axis.y * num3 + axis.z * num4);
    return num * num5;
}

/**
 * Calculate a position between the points specified by current and target, moving 
 * no farther than the distance specified by maxDistanceDelta.
 * @param current The position to move from.
 * @param target The position to move towards.
 * @param maxDistanceDelta Distance to move current per call.
 * @returns The new position.
 */
export function moveTowards(current: Vector3, target: Vector3, maxDistanceDelta: number, out = new Vector3()) {

    const num  = target.x - current.x;
    const num2 = target.y - current.y;
    const num3 = target.z - current.z;
    const num4 = num * num + num2 * num2 + num3 * num3;

    if (num4 == 0 || (maxDistanceDelta >= 0 && num4 <= maxDistanceDelta * maxDistanceDelta)) {
        return assign(out, target);
    }

    const num5 = Math.sqrt(num4);

    return set(
        out,
        current.x + num / num5 * maxDistanceDelta,
        current.y + num2 / num5 * maxDistanceDelta,
        current.z + num3 / num5 * maxDistanceDelta
    );
}

export const overSqrt2 = 0.7071067811865475244008443621048490;

export function orthoNormalVector(vec: Vector3, out = new Vector3()) {

	if (Math.abs(vec.z) > overSqrt2) {			
		const a = vec.y * vec.y + vec.z * vec.z;
		const k = 1 / Math.sqrt(a);
        return set(out, 0, -vec.z * k, vec.y * k);
    }

    const a = vec.x * vec.x + vec.y * vec.y;
    const k = 1 / Math.sqrt(a);

    return set(out, -vec.y * k, vec.x * k, 0);
}

export function clampedMove(lhs: number, rhs: number, clampedDelta: number) {

	const delta = rhs - lhs
	
	if (delta > 0) {
		return lhs + Math.min(delta, clampedDelta);
    }

    return lhs - Math.min(-delta, clampedDelta);
}

const angleAxisTmpV3 = new Vector3();

export function angleAxis(angle: number, axis: Vector3, out = new Quaternion()) {

    const normalizeAxis = normalize(axis, angleAxisTmpV3);
    const halfAngle = angle * Math.deg2Rad * 0.5;
    const s = Math.sin(halfAngle);
    const c = Math.cos(halfAngle);
    
    return QuaternionExBase.set(
        out,
        normalizeAxis.x * s,
        normalizeAxis.y * s,
        normalizeAxis.z * s,
        c
    );
}

/**
 * Rotates a vector current towards target.
 * @param current The vector being managed.
 * @param target The vector.
 * @param maxRadiansDelta The maximum angle in radians allowed for this rotation.
 * @param maxMagnitudeDelta The maximum allowed change in vector magnitude for this rotation.
 * @param out The result vector
 * @returns The location that RotateTowards generates.
 */
export function rotateTowards(current: Vector3, target: Vector3, maxRadiansDelta: number, maxMagnitudeDelta: number, out = new Vector3()) {

	const len1 = magnitude(current);
	const len2 = magnitude(target);
	
	if (len1 > 1e-6 && len2 > 1e-6) {

		const from  = divideScalar(current, len1);
		const to    = divideScalar(target, len2);
		const cosom = dot(from, to);
        
		if (cosom > 1 - 1e-6) {
			return moveTowards(current, target, maxMagnitudeDelta, out);
        }

        if (cosom < -1 + 1e-6) {

            const axis  = orthoNormalVector(from);
			const q     = angleAxis(maxRadiansDelta * Math.rad2Deg, axis);
			const delta = clampedMove(len1, len2, maxMagnitudeDelta);

			       multiplyQV3(q, from, out);
			return multiplyScalarRef(out, delta);
        }

        const angle = Math.acos(cosom);
        const axis  = cross(from, to);

        normalizeRef(axis);

        const q     = angleAxis(Math.min(maxRadiansDelta, angle) * Math.rad2Deg, axis);
        const delta = clampedMove(len1, len2, maxMagnitudeDelta);

               multiplyQV3(q, from, out);
        return multiplyScalarRef(out, delta);
    }

	return moveTowards(current, target, maxMagnitudeDelta, out);
}

export function smoothDamp(current: Vector3, target: Vector3, refCurrentVelocity: Vector3, smoothTime: number, maxSpeed: number, deltaTime: number, out = new Vector3()) {

    let num = 0;
    let num2 = 0;
    let num3 = 0;

    smoothTime = Math.max(0.0001, smoothTime);
    const num4 = 2 / smoothTime;
    const num5 = num4 * deltaTime;
    const num6 = 1 / (1 + num5 + 0.48 * num5 * num5 + 0.235 * num5 * num5 * num5);
    let num7 = current.x - target.x;
    let num8 = current.y - target.y;
    let num9 = current.z - target.z;

    const vector = clone(target);

    const num10 = maxSpeed * smoothTime;
    const num11 = num10 * num10;
    const num12 = num7 * num7 + num8 * num8 + num9 * num9;
    
    if (num12 > num11) {
        const num13 = Math.sqrt(num12);
        num7 = num7 / num13 * num10;
        num8 = num8 / num13 * num10;
        num9 = num9 / num13 * num10;
    }

    target.x = current.x - num7;
    target.y = current.y - num8;
    target.z = current.z - num9;

    const num14 = (refCurrentVelocity.x + num4 * num7) * deltaTime;
    const num15 = (refCurrentVelocity.y + num4 * num8) * deltaTime;
    const num16 = (refCurrentVelocity.z + num4 * num9) * deltaTime;

    refCurrentVelocity.x = (refCurrentVelocity.x - num4 * num14) * num6;
    refCurrentVelocity.y = (refCurrentVelocity.y - num4 * num15) * num6;
    refCurrentVelocity.z = (refCurrentVelocity.z - num4 * num16) * num6;

    num = target.x + (num7 + num14) * num6;
    num2 = target.y + (num8 + num15) * num6;
    num3 = target.z + (num9 + num16) * num6;

    const num17 = vector.x - current.x;
    const num18 = vector.y - current.y;
    const num19 = vector.z - current.z;
    const num20 = num - vector.x;
    const num21 = num2 - vector.y;
    const num22 = num3 - vector.z;

    if (num17 * num20 + num18 * num21 + num19 * num22 > 0) {
        num = vector.x;
        num2 = vector.y;
        num3 = vector.z;
        refCurrentVelocity.x = (num - vector.x) / deltaTime;
        refCurrentVelocity.y = (num2 - vector.y) / deltaTime;
        refCurrentVelocity.z = (num3 - vector.z) / deltaTime;
    }

    return set(out, num, num2, num3);
}

export const Vector3Ex = {
    set,
    assign,
    clone,
    normalize,
    cross,
    reflect,
    clampMagnitude,
    setMagnitude,
    distance,
    distanceSquared,
    dot,
    min,
    max,
    abs,
    clamp,
    lerp,
    lerpUnclamped,
    slerp,
    slerpUnclamped,
    equals,
    notEquals,
    squareRoot,
    add,
    subtract,
    multiply,
    multiplyScalar,
    multiplyScalarL,
    multiplyQV3,
    divide,
    divideScalar,
    divideScalarL,
    negate,
    orthoNormalize,
    angleAxis,
    angle,
    signedAngle,
    sqrMagnitude,
    magnitude,
    project,
    projectOnPlane,
    scale,
    moveTowards,
    rotateTowards,
    smoothDamp,
}

export const Vector3ExByRef = {
    normalizeRef,
    crossRef,
    crossRef2,
    reflectRef,
    minRef,
    maxRef,
    absRef,
    clampRef,
    lerpRef,
    lerpRef2,
    lerpUnclampedRef,
    lerpUnclampedRef2,
    slerpRef,
    slerpRef2,
    slerpUnclampedRef,
    slerpUnclampedRef2,
    addRef,
    addRef2,
    subtractRef,
    subtractRef2,
    multiplyRef,
    multiplyRef2,
    multiplyScalarRef,
    multiplyScalarLRef2,
    multiplyQV3Ref2,
    divideRef,
    divideRef2,
    divideScalarRef,
    divideScalarLRef2,
    negateRef,
    projectRef,
    projectRef2,
    projectOnPlaneRef,
    projectOnPlaneRef2,
    scaleRef,
    scaleRef2,
    clampMagnitudeRef,
    setMagnitudeRef,
}

export const Vec3Ex = {
    ...UnityVector3,
    ...Vector3Ex,
    ...Vector3ExByRef,
};

export default Vec3Ex;
