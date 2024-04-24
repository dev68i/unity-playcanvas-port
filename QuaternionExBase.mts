import { Quaternion } from "./Types.mjs";

/**
 * Sets the specified quaternion to the supplied numerical values.
 *
 * @param x - The x component of the quaternion.
 * @param y - The y component of the quaternion.
 * @param z - The z component of the quaternion.
 * @param w - The w component of the quaternion.
 * @returns Self for chaining.
 * @example
 * const q = new Quaternion();
 * Set(q, 1, 0, 0, 0);
 *
 * // Outputs 1, 0, 0, 0
 * console.log("The result of the quaternion set is: " + q.toString());
 */
export function set(to: Quaternion, x: number, y: number, z: number, w: number) {
    to.x = x;
    to.y = y;
    to.z = z;
    to.w = w;
    return to;
}

/**
 * Copies the contents of a source quaternion to a destination quaternion.
 *
 * @param to - The quaternion to set
 * @param from - A quaternion to copy to the specified quaternion.
 * @returns Self for chaining.
 * @example
 * const src = new Quaternion(10, 20, 30, 1);
 * const dst = new Quaternion();
 *
 * assign(dst, src);
 *
 * console.log("The two quaternions are " + (Equals(src, dst) ? "equal" : "different"));
 */
export function assign(to: Quaternion, from: Readonly<Quaternion>) {
    to.x = from.x;
    to.y = from.y;
    to.z = from.z;
    to.w = from.w;
    return to;
}

/**
 * Returns an identical copy of the specified quaternion.
 *
 * @returns A quaternion containing the result of the cloning.
 * @example
 * const v = new quaternion(10, 20, 30, 1);
 * const vclone = Clone(v);
 * console.log("The result of the cloning is: " + vclone.toString());
 */
export function clone(from: Readonly<Quaternion>) {
    return new Quaternion(from.x, from.y, from.z, from.w);
}

export default {
    set,
    assign,
    clone,
}
