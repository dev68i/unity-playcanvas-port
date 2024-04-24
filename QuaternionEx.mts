/**
 * Helper methods for dealing with Quaternions.
 */
import { Quaternion, Vector3, float } from "./Types.mjs";
import { set, assign, clone } from "./QuaternionExBase.mjs";
import Math, { kEpsilon } from './Math.mjs';
import Vector3Ex, { Vec3Ex } from "./Vector3Ex.mjs";
import DVector3 from "./Vector3Dfl.mjs";

const halfDegToRad = 0.5 * Math.deg2Rad;

/**
 * Returns a boolean indicating whether the two given quaternions are equal.
 */
export function equals(left: Readonly<Quaternion>, right: Readonly<Quaternion>) {
    return  left.x === right.x &&
            left.y === right.y &&
            left.z === right.z &&
            left.w === right.w;
}

/**
 * Returns a boolean indicating whether the two given quaternions are not equal.
 */
export function notEquals(left: Readonly<Quaternion>, right: Readonly<Quaternion>) {
    return  left.x !== right.x ||
            left.y !== right.y ||
            left.z !== right.z ||
            left.w !== right.w;
}

const lookRotationTmp1V3 = new Vector3();
const lookRotationTmp2V3 = new Vector3();
const lookRotationTmp3V3 = new Vector3();

export function lookRotation(forward: Readonly<Vector3>, up: Readonly<Vector3> = DVector3.up, out = new Quaternion()) {
    
    const vector  = Vector3Ex.normalize(forward, lookRotationTmp1V3);
    const vector2 = Vector3Ex.normalizeRef(Vector3Ex.cross(up, vector, lookRotationTmp2V3));
    const vector3 = Vector3Ex.cross(vector, vector2, lookRotationTmp3V3);

    const m00 = vector2.x;
    const m01 = vector2.y;
    const m02 = vector2.z;
    const m10 = vector3.x;
    const m11 = vector3.y;
    const m12 = vector3.z;
    const m20 = vector.x;
    const m21 = vector.y;
    const m22 = vector.z;

    const num8 = (m00 + m11) + m22;

    if (num8 > 0) {

        let num = Math.sqrt(num8 + 1);
        out.w = num * 0.5;
        num = 0.5 / num;
        out.x = (m12 - m21) * num;
        out.y = (m20 - m02) * num;
        out.z = (m01 - m10) * num;

        return out;
    }

    if ((m00 >= m11) && (m00 >= m22)) {

        const num7 = Math.sqrt(((1 + m00) - m11) - m22);
        const num4 = 0.5 / num7;

        out.x = 0.5 * num7;
        out.y = (m01 + m10) * num4;
        out.z = (m02 + m20) * num4;
        out.w = (m12 - m21) * num4;

        return out;
    }

    if (m11 > m22) {

        const num6 = Math.sqrt(((1 + m11) - m00) - m22);
        const num3 = 0.5 / num6;

        out.x = (m10+ m01) * num3;
        out.y = 0.5 * num6;
        out.z = (m21 + m12) * num3;
        out.w = (m20 - m02) * num3;

        return out; 
    }

    const num5 = Math.sqrt(((1 + m22) - m00) - m11);
    const num2 = 0.5 / num5;

    out.x = (m20 + m02) * num2;
    out.y = (m21 + m12) * num2;
    out.z = 0.5 * num5;
    out.w = (m01 - m10) * num2;

    return out;
}

export const multiplyV3 = Vector3Ex.multiplyQV3;
export const multiplyV3Ref2 = Vector3Ex.multiplyQV3Ref2;

export interface IM3 {
    x: number,
    y: number,
    z: number
}

export interface IMatrixOfObj {
    x: Readonly<IM3>,
    y: Readonly<IM3>,
    z: Readonly<IM3>,
}

type TMK = 'x' | 'y' | 'z';
type TMK3 = {x: TMK, y: TMK, z: TMK};

const _next: TMK3 = {x: 'y', y: 'z', z: 'x'};
const _q = new Vector3();

export function matrixToQuaternion(rot: Readonly<IMatrixOfObj>, out = new Quaternion()) {

    const trace = rot.x.x + rot.y.y + rot.z.z;

	if (trace > 0) {

		let s = Math.sqrt(trace + 1);
		out.w = 0.5 * s;
		s = 0.5 / s;
		out.x = (rot.z.y - rot.y.z) * s;
		out.y = (rot.x.z - rot.z.x) * s;
		out.z = (rot.y.x - rot.x.y) * s;
        return normalizeRef(out);
    }

    let i: TMK = 'x';

    if (rot.y.y > rot.x.x) {
        i = 'y';
    }

    if (rot.z.z > rot[i][i]) {
        i = 'z';
    }

    const j = _next[i];
    const k = _next[j];
		
    const t = rot[i][i] - rot[j][j] - rot[k][k] + 1;
    const s = 0.5 / Math.sqrt(t);

    _q[i] = s * t;

    const w = (rot[k][j] - rot[j][k]) * s;
      _q[j] = (rot[j][i] + rot[i][j]) * s;
      _q[k] = (rot[k][i] + rot[i][k]) * s;
		
    set(out, _q.x, _q.y, _q.z, w);

    return normalizeRef(out);
}

// Test not success
/*
export function fromToRotation(from: Vector3, to: Vector3) {
    const axis  = Cross(from, to);
    const angle = Angle(from, to);
    return AngleAxis(angle, axis.normalize());
}
*/

const fromToRotationTmp1V3 = new Vector3();
const fromToRotationTmp2V3 = new Vector3();
const fromToRotationTmp3V3 = new Vector3();
const fromToRotationTmp4V3 = new Vector3();

const fromToRotationV3Tmp2Obj = { x: 0, y: 0, z: 0};
const fromToRotationV3Tmp3Obj = { x: 0, y: 0, z: 0};
const fromToRotationV3Tmp4Obj = { x: 0, y: 0, z: 0};
const fromToRotationV3TmpMatrix = {
    x: fromToRotationV3Tmp2Obj,
    y: fromToRotationV3Tmp3Obj,
    z: fromToRotationV3Tmp4Obj
};

export function fromToRotationV3(from: Readonly<Vector3>, to: Readonly<Vector3>, out = new Quaternion()) {

    const nFrom = Vector3Ex.normalize(from, fromToRotationTmp1V3);
	const nTo   = Vector3Ex.normalize(to, fromToRotationTmp2V3);
	
	const e = Vector3Ex.dot(nFrom, nTo);
	
	if (e > 1 - 1e-6) {
		return set(out, 0, 0, 0, 1);
    }

	if (e < -1 + 1e-6) {

		let left = Vec3Ex.set(fromToRotationTmp3V3, 0, nFrom.z, nFrom.y);
		let mag = left.y * left.y + left.z * left.z; //+ left[0] * left[0] = 0
		
		if (mag < 1e-6) {
			left.x = -nFrom.z
			left.y = 0
			left.z = nFrom.x
			mag = left.x * left.x + left.z * left.z;
        }
				
		const invlen = 1 / Math.sqrt(mag);

		left.x *= invlen;
		left.y *= invlen;
		left.z *= invlen;
		
		const up = Vec3Ex.set(fromToRotationTmp4V3, 0, 0, 0);

		up.x = left.y * nFrom.z - left.z * nFrom.y;
		up.y = left.z * nFrom.x - left.x * nFrom.z;
		up.z = left.x * nFrom.y - left.y * nFrom.x;

		const fxx = -nFrom.x * nFrom.x;
		const fyy = -nFrom.y * nFrom.y;
		const fzz = -nFrom.z * nFrom.z;
		
		const fxy = -nFrom.x * nFrom.y;
		const fxz = -nFrom.x * nFrom.z;
		const fyz = -nFrom.y * nFrom.z;

		const uxx = up.x * up.x;
		const uyy = up.y * up.y;
		const uzz = up.z * up.z;
		const uxy = up.x * up.y;
		const uxz = up.x * up.z;
		const uyz = up.y * up.z;

		const lxx = -left.x * left.x;
		const lyy = -left.y * left.y;
		const lzz = -left.z * left.z;
		const lxy = -left.x * left.y;
		const lxz = -left.x * left.z;
		const lyz = -left.y * left.z;
		
        /*
		const rot = 
		[
			[fxx + uxx + lxx, fxy + uxy + lxy, fxz + uxz + lxz],
			[fxy + uxy + lxy, fyy + uyy + lyy, fyz + uyz + lyz],
			[fxz + uxz + lxz, fyz + uyz + lyz, fzz + uzz + lzz],
        ];
        */
        fromToRotationV3Tmp2Obj.x = fxx + uxx + lxx;
        fromToRotationV3Tmp2Obj.y = fxy + uxy + lxy;
        fromToRotationV3Tmp2Obj.z = fxz + uxz + lxz;
        fromToRotationV3Tmp3Obj.x = fxy + uxy + lxy;
        fromToRotationV3Tmp3Obj.y = fyy + uyy + lyy;
        fromToRotationV3Tmp3Obj.z = fyz + uyz + lyz;
        fromToRotationV3Tmp4Obj.x = fxz + uxz + lxz;
        fromToRotationV3Tmp4Obj.y = fyz + uyz + lyz;
        fromToRotationV3Tmp4Obj.z = fzz + uzz + lzz;

		return matrixToQuaternion(fromToRotationV3TmpMatrix, out);
	}

    const v = Vector3Ex.crossRef(nFrom, nTo);
    const h = (1 - e) / Vector3Ex.dot(v, v);
    
    const hx = h * v.x
    const hz = h * v.z
    const hxy = hx * v.y
    const hxz = hx * v.z
    const hyz = hz * v.y
    
    /*
    const rot = 
    [ 					
        [e + hx*v.x, 	hxy - v.z, 		hxz + v.y],
        [hxy + v.z,  	e + h*v.y*v.y, 	hyz-v.x],
        [hxz - v.y,  	hyz + v.x,    	e + hz*v.z],
    ]
    */

    fromToRotationV3Tmp2Obj.x = e + hx*v.x;
    fromToRotationV3Tmp2Obj.y = hxy - v.z;
    fromToRotationV3Tmp2Obj.z = hxz + v.y;
    fromToRotationV3Tmp3Obj.x = hxy + v.z;
    fromToRotationV3Tmp3Obj.y = e + h*v.y*v.y;
    fromToRotationV3Tmp3Obj.z = hyz-v.x;
    fromToRotationV3Tmp4Obj.x = hxz - v.y;
    fromToRotationV3Tmp4Obj.y = hyz + v.x;
    fromToRotationV3Tmp4Obj.z = e + hz*v.z;
    
    return matrixToQuaternion(fromToRotationV3TmpMatrix, out);
}

/**
 * Multiplies two Quaternions together.
 */
export function multiply(value1: Readonly<Quaternion>, value2: Readonly<Quaternion>, out = new Quaternion()) {

    const q1x = value1.x;
    const q1y = value1.y;
    const q1z = value1.z;
    const q1w = value1.w;

    const q2x = value2.x;
    const q2y = value2.y;
    const q2z = value2.z;
    const q2w = value2.w;

    // cross(av, bv)
    const cx = q1y * q2z - q1z * q2y;
    const cy = q1z * q2x - q1x * q2z;
    const cz = q1x * q2y - q1y * q2x;

    const dot = q1x * q2x + q1y * q2y + q1z * q2z;

    return set(
        out,
        q1x * q2w + q2x * q1w + cx,
        q1y * q2w + q2y * q1w + cy,
        q1z * q2w + q2z * q1w + cz,
        q1w * q2w - dot
    );
}

/**
 * Multiplies two Quaternions together.
 */
export function multiplyRef(refValue1: Quaternion, value2: Readonly<Quaternion>) {
    return multiply(refValue1, value2, refValue1);
}

/**
 * Multiplies two Quaternions together.
 */
export function multiplyRef2(value1: Readonly<Quaternion>, refValue2: Quaternion) {
    return multiply(value1, refValue2, refValue2);
}

/**
 * Multiplies a Quaternion by a scalar value.
 */
export function multiplyScalar(value1: Readonly<Quaternion>, value2: number, out = new Quaternion()) {
    return set(
        out,
        value1.x * value2,
        value1.y * value2,
        value1.z * value2,
        value1.w * value2
    );
}

/**
 * Multiplies a Quaternion by a scalar value.
 */
export function multiplyScalarRef(refValue1: Quaternion, value2: number) {
    return multiplyScalar(refValue1, value2, refValue1);
}

export function normalize(value: Readonly<Quaternion>, out = new Quaternion()) {

    const ls = value.x * value.x + value.y * value.y + value.z * value.z + value.w * value.w;
    
    if(ls !== 0) {
        const invNorm = 1.0 / Math.sqrt(ls);
        return set(
            out,
            value.x * invNorm,
            value.y * invNorm,
            value.z * invNorm,
            value.w * invNorm
        );
    }

    return assign(out, value);
}

export function normalizeRef(refValue: Quaternion) {
    return normalize(refValue, refValue);
}

export function inverse(value: Readonly<Quaternion>, out = new Quaternion()) {

    //  -1   (       a              -v       )
    // q   = ( -------------   ------------- )
    //       (  a^2 + |v|^2  ,  a^2 + |v|^2  )
    
    const ls: float = value.x * value.x + value.y * value.y + value.z * value.z + value.w * value.w;

    if(ls !== 0) {
        
        const invNorm: float = 1.0 / ls;

        return set(
            out,
            -value.x * invNorm,
            -value.y * invNorm,
            -value.z * invNorm,
            value.w * invNorm
        );
    }

    return assign(out, value);
}

export function inverseRef(refValue: Quaternion) {
    return inverse(refValue, refValue);
}

export function lerpUnclamped(quaternion1: Readonly<Quaternion>, quaternion2: Readonly<Quaternion>, amount: number, out = new Quaternion()) {

    const t = amount;
    const t1 = 1.0 - t;
    const q1x = quaternion1.x;
    const q1y = quaternion1.y;
    const q1z = quaternion1.z;
    const q1w = quaternion1.w;
    const q2x = quaternion2.x;
    const q2y = quaternion2.y;
    const q2z = quaternion2.z;
    const q2w = quaternion2.w;

    const dot = q1x * q2x + q1y * q2y +
                q1z * q2z + q1w * q2w;
    
    if (dot >= 0)
    {
        out.x = t1 * q1x + t * q2x;
        out.y = t1 * q1y + t * q2y;
        out.z = t1 * q1z + t * q2z;
        out.w = t1 * q1w + t * q2w;
    }
    else
    {
        out.x = t1 * q1x - t * q2x;
        out.y = t1 * q1y - t * q2y;
        out.z = t1 * q1z - t * q2z;
        out.w = t1 * q1w - t * q2w;
    }

    return normalizeRef(out);
}

export function slerpUnclamped(quaternion1: Readonly<Quaternion>, quaternion2: Readonly<Quaternion>, amount: number, out = new Quaternion()) {

    // Algorithm sourced from:
    // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
    const lx = quaternion1.x;
    const ly = quaternion1.y;
    const lz = quaternion1.z;
    const lw = quaternion1.w;

    let rx = quaternion2.x;
    let ry = quaternion2.y;
    let rz = quaternion2.z;
    let rw = quaternion2.w;

    // Calculate angle between them.
    let cosHalfTheta = lw * rw + lx * rx + ly * ry + lz * rz;

    if (cosHalfTheta < 0) {
        rw = -rw;
        rx = -rx;
        ry = -ry;
        rz = -rz;
        cosHalfTheta = -cosHalfTheta;
    }

    // If lhs == rhs or lhs == -rhs then theta == 0 and we can return lhs
    if (Math.abs(cosHalfTheta) >= 1) {
        return set(out, lx, ly, lz, lw);
    }

    // Calculate temporary values.
    const halfTheta    = Math.acos(cosHalfTheta);
    const sinHalfTheta = Math.sqrt(1 - cosHalfTheta * cosHalfTheta);

    // If theta = 180 degrees then result is not fully defined
    // we could rotate around any axis normal to qa or qb
    if (Math.abs(sinHalfTheta) < 0.001) {
        return set(
            out,
            (lx * 0.5 + rx * 0.5),
            (ly * 0.5 + ry * 0.5),
            (lz * 0.5 + rz * 0.5),
            (lw * 0.5 + rw * 0.5)
        );
    }

    const ratioA = Math.sin((1 - amount) * halfTheta) / sinHalfTheta;
    const ratioB = Math.sin(amount * halfTheta) / sinHalfTheta;
    
    // Calculate Quaternion.
    return set(
        out,
        (lx * ratioA + rx * ratioB),
        (ly * ratioA + ry * ratioB),
        (lz * ratioA + rz * ratioB),
        (lw * ratioA + rw * ratioB)
    );
}

/**
 * Interpolates between fromRotation and toRotation by weight and normalizes the result afterwards. The parameter t is clamped to the range [0, 1].
 * @param fromRotation Start value, returned when waight = 0.
 * @param toRotation End value, returned when weight = 1.
 * @param weight Interpolation ratio.
 */
export function lerp(fromRotation: Readonly<Quaternion>, toRotation: Readonly<Quaternion>, weight: number, out = new Quaternion()) {
    if (weight <= 0) return assign(out, fromRotation);
    if (weight >= 1) return assign(out, toRotation);
    return lerpUnclamped(fromRotation, toRotation, weight, out);
}

/**
 * Optimized Quaternion.lerp
 */
export const optLerp = lerp;

/**
 * Interpolates between fromRotation and toRotation by weight and normalizes the result afterwards. The parameter t is clamped to the range [0, 1].
 * @param refFromRotation Ref of start value, returned when waight = 0.
 * @param toRotation End value, returned when weight = 1.
 * @param weight Interpolation ratio.
 */
export function lerpRef(refFromRotation: Quaternion, toRotation: Readonly<Quaternion>, weight: number) {
    return lerp(refFromRotation, toRotation, weight, refFromRotation);
}

/**
 * Interpolates between fromRotation and toRotation by weight and normalizes the result afterwards. The parameter t is clamped to the range [0, 1].
 * @param fromRotation Start value, returned when waight = 0.
 * @param toRotation End value, returned when weight = 1.
 * @param weight Interpolation ratio.
 */
export function lerpRef2(fromRotation: Readonly<Quaternion>, refToRotation: Quaternion, weight: number) {
    return lerp(fromRotation, refToRotation, weight, refToRotation);
}

/**
 * Spherically interpolates between quaternions fromRotation and toRotation by ratio weight. The parameter weight is clamped to the range [0, 1].
 * @param fromRotation Start value, returned when waight = 0.
 * @param toRotation End value, returned when weight = 1.
 * @param weight Interpolation ratio.
 */
export function slerp(fromRotation: Readonly<Quaternion>, toRotation: Readonly<Quaternion>, weight: number, out = new Quaternion()) {
    if (weight <= 0) return assign(out, fromRotation);
    if (weight >= 1) return assign(out, toRotation);
    return slerpUnclamped(fromRotation, toRotation, weight, out);
}

/**
 * Spherically interpolates between quaternions fromRotation and toRotation by ratio weight. The parameter weight is clamped to the range [0, 1].
 * @param fromRotation Start value, returned when waight = 0.
 * @param toRotation End value, returned when weight = 1.
 * @param weight Interpolation ratio.
 */
export function slerpRef(refFromRotation: Quaternion, toRotation: Readonly<Quaternion>, weight: number) {
    return slerp(refFromRotation, toRotation, weight, refFromRotation);
}

/**
 * Spherically interpolates between quaternions fromRotation and toRotation by ratio weight. The parameter weight is clamped to the range [0, 1].
 * @param fromRotation Start value, returned when waight = 0.
 * @param toRotation End value, returned when weight = 1.
 * @param weight Interpolation ratio.
 */
export function slerpRef2(fromRotation: Readonly<Quaternion>, refRoRotation: Quaternion, weight: number) {
    return slerp(fromRotation, refRoRotation, weight, refRoRotation);
}

// axis must be unit, rad is angle in radians
export function axisAngle(axis: Readonly<Vector3>, rad: number, out = new Quaternion()) {
    rad *= .5;
    const scalar = Math.sin(rad);
    return set(
        out,
        axis.x * scalar,
        axis.y * scalar,
        axis.z * scalar,
        Math.cos(rad)
    );
}

export function rotateTowards(from: Readonly<Quaternion>, to: Readonly<Quaternion>, maxDegreesDelta: number, out = new Quaternion()) {
	const angleValue = angle(from, to)
	if (angleValue == 0) {
		return assign(out, to);
    }
	const t = Math.min(1, maxDegreesDelta / angleValue);
	return slerpUnclamped(from, to, t, out);
}

/**
 * Returns yaw angle (-180 - 180) of 'forward' vector relative to rotation space defined by spaceForward and spaceUp axes.
 */
export function getYawV3(space: Readonly<Quaternion>, forward: Readonly<Vector3>) {
    const dirLocal = multiplyV3(inverse(space), forward);
    return Math.atan2(dirLocal.x, dirLocal.z) * Math.rad2Deg;
}

/**
 * Returns pitch angle (-90 - 90) of 'forward' vector relative to rotation space defined by spaceForward and spaceUp axes.
 */
export function getPitchV3(space: Readonly<Quaternion>, forward: Readonly<Vector3>) {
    const normalizeForward = Vector3Ex.normalize(forward);
    const dirLocal = multiplyV3(inverse(space), normalizeForward);
    return -Math.asin(dirLocal.y) * Math.rad2Deg;
}

/**
 * Returns bank angle (-180 - 180) of 'forward' and 'up' vectors relative to rotation space defined by spaceForward and spaceUp axes.
 */
export function getBankV3(space: Readonly<Quaternion>, forward: Readonly<Vector3>, up: Vector3) {
    const spaceUp    = multiplyV3(space, DVector3.up);
    const invSpace   = inverse(space);
    const newForward = multiplyV3(invSpace, forward);
    const newUp      = multiplyV3(invSpace, up);
    const q          = inverseRef(lookRotation(spaceUp, newForward));
    const newUp2     = multiplyV3(q, newUp);
    return Math.atan2(newUp2.x, newUp2.z) * Math.rad2Deg;
}

/**
 * Returns yaw angle (-180 - 180) of 'forward' vector relative to rotation space defined by spaceForward and spaceUp axes.
 */
export function getYaw(space: Readonly<Quaternion>, rotation: Readonly<Quaternion>) {
    const dirLocal = multiplyV3Ref2(inverse(space), multiplyV3(rotation, DVector3.forward));
    return Math.atan2(dirLocal.x, dirLocal.z) * Math.rad2Deg;
}

/**
 * Returns pitch angle (-90 - 90) of 'forward' vector relative to rotation space defined by spaceForward and spaceUp axes.
 */
export function getPitch(space: Readonly<Quaternion>, rotation: Readonly<Quaternion>) {
    const dirLocal = multiplyV3Ref2(inverse(space), multiplyV3(rotation, DVector3.forward));
    return -Math.asin(dirLocal.y) * Math.rad2Deg;
}

/*
 * Returns bank angle (-180 - 180) of 'forward' and 'up' vectors relative to rotation space defined by spaceForward and spaceUp axes.
 */
export function getBank(space: Readonly<Quaternion>, rotation: Readonly<Quaternion>) {
    const spaceUp  = multiplyV3(space, DVector3.up);
    const invSpace = inverse(space);
    const forward  = multiplyV3(invSpace, multiplyV3(rotation, DVector3.forward));
    const up       = multiplyV3(invSpace, multiplyV3(rotation, DVector3.up));
    const q        = inverseRef(lookRotation(spaceUp, forward));
    const newUp    = multiplyV3(q, up);
    return Math.atan2(newUp.x, newUp.z) * Math.rad2Deg;
}

/**
 * Returns the rotation from identity Quaternion to "q", interpolated linearily by "weight".
 */
export function linearBlend(q: Readonly<Quaternion>, weight: number, out = new Quaternion()) {
    if (weight <= 0) return assign(out, Quaternion.IDENTITY);
    if (weight >= 1) return assign(out, q);
    return lerp(Quaternion.IDENTITY, q, weight, out);
}

/**
 * Returns the rotation from identity Quaternion to "q", interpolated spherically by "weight".
 */
export function sphericalBlend(q: Readonly<Quaternion>, weight: number, out = new Quaternion()) {
    if (weight <= 0) return assign(out, Quaternion.IDENTITY);
    if (weight >= 1) return assign(out, q);
    return slerp(Quaternion.IDENTITY, q, weight, out);
}

export function toAxisAngleRad(quaternion: Readonly<Quaternion>, out = new Vector3()): [number, Vector3] {

    let q = Math.abs(quaternion.w) > 1.0
        ? normalize(quaternion)
        : quaternion;

    const angle = 2.0 * Math.acos(q.w); // angle
    const den = Math.sqrt(1.0 - q.w * q.w);

    const axis = den > 0.0001
        ? Vector3Ex.set(out, q.x / den, q.y / den, q.z / den)
        : Vector3Ex.set(out, 1, 0, 0);
        // ^
        // This occurs when the angle is zero.
        // Not a problem: just set an arbitrary normalized axis.

    return [
        angle,
        axis
    ]
}

export function toAngleAxis(q: Readonly<Quaternion>, out = new Vector3()): [number, Vector3] {
    const [angle, axis] = toAxisAngleRad(q, out);
    return [
        angle * Math.rad2Deg,
        axis,
    ]
}

/**
 * Creates a FromToRotation, but makes sure its axis remains fixed near to the Quaternion singularity point.
 * @param fromDirection From direction.
 * @param toDirection To direction.
 * @param axis Axis. Should be normalized before passing into this method.
 * @returns The from to rotation around an axis.
 */
export function fromToAroundAxis(fromDirection: Readonly<Vector3>, toDirection: Readonly<Vector3>, axis: Readonly<Vector3>, out = new Quaternion()) {

    const fromTo = fromToRotationV3(fromDirection, toDirection, out);
    let [ angle = 0, freeAxis = DVector3.zero ] = toAngleAxis(fromTo);
    
    const dot = Vector3Ex.dot(freeAxis, axis);

    if (dot < 0) {
        angle = -angle;
    }
    
    return Vector3Ex.angleAxis(angle, axis, out);
}

const rotationTmp1Qt = new Quaternion();

/**
 * Gets the rotation that can be used to convert a rotation from one axis space to another.
 */
export function rotationToLocalSpace(space: Readonly<Quaternion>, rotation: Readonly<Quaternion>, out = new Quaternion()) {
    const t = inverse(space, rotationTmp1Qt);
    multiplyRef(t, rotation);
    return inverse(t, out);
}

/**
 * Gets the Quaternion from rotation "from" to rotation "to".
 */
export function fromToRotation(from: Readonly<Quaternion>, to: Readonly<Quaternion>, out = new Quaternion()) {

    if (equals(to, from)) {
        return assign(out, Quaternion.IDENTITY);
    }

    return multiply(to, inverse(from, rotationTmp1Qt), out);
}

export function isEqualUsingDot(dot: number) {
    // Returns false in the presence of NaN values.
    return dot > 1.0 - kEpsilon;
}

export function dot(a: Readonly<Quaternion>, b: Readonly<Quaternion>) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

export function angle(a: Readonly<Quaternion>, b: Readonly<Quaternion>) {
    const dotValue = Math.min(Math.abs(dot(a, b)), 1.0);
    return isEqualUsingDot(dotValue) ? 0.0 : Math.acos(dotValue) * 2.0 * Math.rad2Deg;
}

const angleAxisTmpV3 = new Vector3();

export function angleAxis(angle: number, axis: Readonly<Vector3>, out = new Quaternion()) {

    const normAxis = Vector3Ex.normalize(axis, angleAxisTmpV3);
    angle = angle * halfDegToRad;
    const s = Math.sin(angle);
    
    const w = Math.cos(angle);
    const x = normAxis.x * s;
    const y = normAxis.y * s;
    const z = normAxis.z * s;
    
    return set(out, x, y, z, w);
}

/**
 * Clamps the rotation similar to V3Tools.clampDirection.
 */
export function clampRotation(rotation: Readonly<Quaternion>, clampWeight: number, clampSmoothing: number, out = new Quaternion()) {

    if (clampWeight >= 1) return assign(out, Quaternion.IDENTITY);
    if (clampWeight <= 0) return assign(out, rotation);

    const angleValue = angle(Quaternion.IDENTITY, rotation);
    const dot = 1 - (angleValue / 180);
    const targetClampMlp = Math.clamp(1 - ((clampWeight - dot) / (1 - dot)), 0, 1);

    let clampMlp = Math.clamp(dot / clampWeight, 0, 1);
    
    // Sine smoothing iterations
    for (let i = 0; i < clampSmoothing; i++) {
        const sinF = clampMlp * Math.PI * 0.5;
        clampMlp = Math.sin(sinF);
    }
    
    return slerp(Quaternion.IDENTITY, rotation, clampMlp * targetClampMlp, out);
}

export function clampAngle(angle: number, clampWeight: number, clampSmoothing: number) {

    if (clampWeight >= 1) return 0;
    if (clampWeight <= 0) return angle;

    const dot = 1 - (Math.abs(angle) / 180);
    const targetClampMlp = Math.clamp(1 - ((clampWeight - dot) / (1 - dot)), 0, 1);

    let clampMlp = Math.clamp(dot / clampWeight, 0, 1);
    
    // Sine smoothing iterations
    for (let i = 0; i < clampSmoothing; i++) {
        const sinF = clampMlp * Math.PI * 0.5;
        clampMlp = Math.sin(sinF);
    }
    
    return Math.lerp(0, angle, clampMlp * targetClampMlp);
}

/**
 * Used for matching the rotations of objects that have different orientations.
 */
export function matchRotation(
    targetRotation: Readonly<Quaternion>,
    targetAxis1: Readonly<Vector3>,
    targetAxis2: Readonly<Vector3>,
    axis1: Readonly<Vector3>,
    axis2: Readonly<Vector3>,
    out = new Quaternion()
) {
    const f = lookRotation(axis1, axis2);
    const t = lookRotation(targetAxis1, targetAxis2);
    multiplyRef(t, targetRotation);
    inverseRef(f);
    return multiply(t, f, out);
}

/**
 * Converts an angular value from 0 to 360 representation to -180 to 180.
 * @returns 
 */
export function toBiPolarA(angle: number) {

    angle %= 360;

    if (angle >= 180)  return angle - 360;
    if (angle <= -180) return angle + 360;

    return angle;
}

/**
 * Converts an Euler rotation from 0 to 360 representation to -180 to 180.
 */
export function toBiPolar(euler: Readonly<Vector3>, out = new Vector3()) {
    return Vector3Ex.set(
        out,
        toBiPolarA(euler.x),
        toBiPolarA(euler.y),
        toBiPolarA(euler.z)
    );
}

/**
 * Mirrors a Quaternion on the YZ plane in provided rotation space.
 * Mirrors a Quaternion on the world space YZ plane.
 */
export function mirrorYZ(r: Readonly<Quaternion>, space?: Readonly<Quaternion>, out = new Quaternion()) {

    if(space) {

        const newR = inverse(space);

        multiplyRef(newR, r);
        
        const forward = multiplyV3(newR, DVector3.forward);
        const up      = multiplyV3(newR, DVector3.up);

        forward.x *= -1;
        up.x *= -1;

        lookRotation(forward, up, newR);

        return multiply(newR, space, out);
    }

    const forward = multiplyV3(r, DVector3.forward);
    const up      = multiplyV3(r, DVector3.up);

    forward.x *= -1;
    up.x *= -1;

    return lookRotation(forward, up, out);
}

/**
 * Returns a rotation that rotates z degrees around the z axis, x degrees around
 */
export function eulerXYZ(x: number, y: number, z: number, out = new Quaternion()): Quaternion {

	x *= halfDegToRad;
    y *= halfDegToRad;
    z *= halfDegToRad;
	
	const sinX = Math.sin(x);
    const cosX = Math.cos(x);
    const sinY = Math.sin(y);
    const cosY = Math.cos(y);
    const sinZ = Math.sin(z);
    const cosZ = Math.cos(z);
        
    return set(
        out,
        cosY * sinX * cosZ + sinY * cosX * sinZ,
        sinY * cosX * cosZ - cosY * sinX * sinZ,
        cosY * cosX * sinZ - sinY * sinX * cosZ,
        cosY * cosX * cosZ + sinY * sinX * sinZ
    );
}

/**
 * Returns a rotation that rotates z degrees around the z axis, x degrees around
 */
export function euler(euler: Readonly<Vector3>, out?: Quaternion): Quaternion {
    return eulerXYZ(euler.x, euler.y, euler.z, out);
}

export function toEulerRadAngles(q: Readonly<Quaternion>, out: Vector3 = new Vector3()): Vector3 {

    const qx = q.x;
    const qy = q.y;
    const qz = q.z;
    const qw = q.w;

    const a2 = 2 * (qw * qy - qx * qz);

    if (a2 <= -0.99999) {
        out.x = 2 * Math.atan2(qx, qw);
        out.y = -Math.PI / 2;
        out.z = 0;
    } else if (a2 >= 0.99999) {
        out.x = 2 * Math.atan2(qx, qw);
        out.y = Math.PI / 2;
        out.z = 0;
    } else {
        out.x = Math.atan2(2 * (qw * qx + qy * qz), 1 - 2 * (qx * qx + qy * qy));
        out.y = Math.asin(a2);
        out.z = Math.atan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy * qy + qz * qz));
    }

    return out;
}

export function toEulerAngles(q: Readonly<Quaternion>, out: Vector3 = new Vector3()): Vector3 {
    toEulerRadAngles(q, out);
    out.x *= Math.rad2Deg;
    out.y *= Math.rad2Deg;
    out.z *= Math.rad2Deg;
    return out;
}

/**
 * Clamping Euler angles
 */
export function clampEulerAngle(angle: number, min: number, max: number) {
    if (angle < -360) angle += 360;
    if (angle > 360) angle -= 360;
    return Math.clamp(angle, min, max);
}

export const QuaternionEx = {
    set,
    assign,
    clone,
    equals,
    notEquals,
    lookRotation,
    inverse,
    multiply,
    multiplyV3,
    fromToRotation,
    fromToRotationV3,
    lerpUnclamped,
    slerpUnclamped,
    lerp,
    slerp,
    optLerp,
    getYaw,
    getPitch,
    getBank,
    getYawV3,
    getPitchV3,
    getBankV3,
    clampRotation,
    clampAngle,
    matchRotation,
    toBiPolarA,
    toBiPolar,
    mirrorYZ,
    euler,
    eulerXYZ,
    clampEulerAngle,
    angle,
    angleAxis,
    axisAngle,
    rotateTowards,
    rotationToLocalSpace,
    dot,
    normalize,
    toEulerAngles,
}

export const QuaternionExRef = {
    normalizeRef,
    multiplyRef,
    multiplyRef2,
    inverseRef,
    multiplyV3Ref2,
    lerpRef,
    lerpRef2,
    slerpRef,
    slerpRef2,
}

export const QuatEx = {
    ...QuaternionEx,
    ...QuaternionExRef,
}

export default QuatEx;
