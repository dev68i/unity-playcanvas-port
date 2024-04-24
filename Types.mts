export type Vector2 = pc.Vec2 & { z?: never; w?: never; };
export type Vector3 = pc.Vec3 & { w?: never; };
export type Quaternion = pc.Quat;
export type Matrix3x3 = pc.Mat3;
export type Matrix4x4 = pc.Mat4;

export const Quaternion = pc.Quat;
export const Vector2 = pc.Vec2;
export const Vector3 = pc.Vec3;
export const Matrix3x3 = pc.Mat3;
export const Matrix4x4 = pc.Mat4;

export type CameraComponent = pc.CameraComponent;
export const CameraComponent = pc.CameraComponent;

export type GraphNode = pc.GraphNode;
export const GraphNode = pc.GraphNode;

export type Entity = pc.Entity;
export const Entity = pc.Entity;

export type Animator = pc.AnimComponent;
export const Animator = pc.AnimComponent;

export type Color = pc.Color;
export const Color = pc.Color;

export type Float = number;
export type float = number;
export type int = number;
export type Int = number;

export type Curve = pc.Curve;
export const Curve = pc.Curve;

export type AnimationCurve = pc.AnimCurve;
export const AnimationCurve = pc.AnimCurve;

export type PhysicsSystem = pc.RigidBodyComponentSystem;
export const PhysicsSystem = pc.RigidBodyComponentSystem;
export type RaycastResult = pc.RaycastResult;
export const RaycastResult = pc.RaycastResult;

export const RigidBodyComponent = pc.RigidBodyComponent;
export type RigidBodyComponent = pc.RigidBodyComponent;

export const CollisionComponent = pc.CollisionComponent;
export type CollisionComponent = pc.CollisionComponent;

export const ScriptComponent = pc.ScriptComponent;
export type ScriptComponent = pc.ScriptComponent;

export interface Ref<T> {
    value: T
}

export interface Out<T> {
    value: T
}

export type RefObject<T extends object> = T;
export type OutObject<T extends object> = T;

export type PropType<TObj, TProp extends keyof TObj> = TObj[TProp];

export type MaybeTransform  = Entity & { enabled?: boolean; };
export const MaybeTransform = Entity;

export declare type Mutable<T extends object> = {
    -readonly [K in keyof T]: T[K]
}

export interface IVector2 {
    x: number,
    y: number,
    z?: never,
    w?: never,
}

export interface IVector3 {
    x: number,
    y: number,
    z: number,
    w?: never,
}

export interface IQuternion {
    x: number,
    y: number,
    z: number,
    w: number,
}
