import { Vector3 } from "./Types.mjs"

export const UnityVector3 = {

    zero:    Object.freeze(new Vector3(0, 0, 0)),
    one:     Object.freeze(new Vector3(1, 1, 1)),
    up:      Object.freeze(new Vector3(0, 1, 0)),
    down:    Object.freeze(new Vector3(0, -1, 0)),
    left:    Object.freeze(new Vector3(-1, 0, 0)),
    right:   Object.freeze(new Vector3(1, 0, 0)),
    forward: Object.freeze(new Vector3(0, 0, 1)),
    back:    Object.freeze(new Vector3(0, 0, -1)),

    positiveInfinity: Object.freeze(new Vector3(globalThis.Number.POSITIVE_INFINITY, globalThis.Number.POSITIVE_INFINITY, globalThis.Number.POSITIVE_INFINITY)),
    negativeInfinity: Object.freeze(new Vector3(globalThis.Number.NEGATIVE_INFINITY, globalThis.Number.NEGATIVE_INFINITY, globalThis.Number.NEGATIVE_INFINITY)),
}

export const PortFromUnityTest = false; // SET TRUE IF UNITY PORTING
export const DVector3 = {
    
    zero:    PortFromUnityTest ? UnityVector3.zero    : Vector3.ZERO,
    one:     PortFromUnityTest ? UnityVector3.one     : Vector3.ONE,
    forward: PortFromUnityTest ? UnityVector3.forward : Vector3.FORWARD,
    back:    PortFromUnityTest ? UnityVector3.back    : Vector3.BACK,
    up:      PortFromUnityTest ? UnityVector3.up      : Vector3.UP,
    down:    PortFromUnityTest ? UnityVector3.down    : Vector3.DOWN,
    left:    PortFromUnityTest ? UnityVector3.left    : Vector3.LEFT,
    right:   PortFromUnityTest ? UnityVector3.right   : Vector3.RIGHT,

    positiveInfinity: UnityVector3.positiveInfinity,
    negativeInfinity: UnityVector3.negativeInfinity
}

export default DVector3;
