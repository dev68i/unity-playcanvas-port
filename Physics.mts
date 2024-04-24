import { CollisionComponent, Matrix4x4, MaybeTransform, Out, OutObject, PhysicsSystem, Quaternion, RaycastResult, RigidBodyComponent, Vector3 } from './Types.mjs';
import Vec3Ex from './Vector3Ex.mjs';
import QuatEx from './QuaternionEx.mjs';
import DVector3, { UnityVector3 } from './Vector3Dfl.mjs';

function ammoHack<T>(callback: () => T) {
    if (!!window['Ammo' as any]) {
        return callback();
    }
    return undefined as T;
}

const helpers = ammoHack(() => ({
    _ammoHalfExtents: new Ammo.btVector3(),
    _ammoRotFrom: new Ammo.btQuaternion(0, 0, 0, 1),
    _ammoRotTo: new Ammo.btQuaternion(0, 0, 0, 1),
    _ammoPosFrom: new Ammo.btVector3(),
    _ammoPosTo: new Ammo.btVector3(),
    _ammoTransformFrom: new Ammo.btTransform(),
    _ammoTransformTo: new Ammo.btTransform(),
    _ammoBlankSphereRadius: 0.5,
    _ammoBlankSphere: new Ammo.btSphereShape(0.5),
    _ammoGhostObjectA: new Ammo.btPairCachingGhostObject(),
    _ammoGhostObjectAWorld: undefined as Ammo.btCollisionWorld | undefined,
    _ammoGhostObjectB: new Ammo.btPairCachingGhostObject(),
    _ammoGhostObjectBWorld: undefined as Ammo.btCollisionWorld | undefined,
    _ammoOverlapContactResponse: new Ammo.ConcreteContactResultCallback(),
    _pcVec31: new Vector3(),
    _pcVec32: new Vector3(),
    _pcVec33: new Vector3(),
    _pcQuat1: new Quaternion(),
    _pcQuat2: new Quaternion(),
    _pcMat4: new Matrix4x4(),
    _pcStartRot: new Quaternion(),
    _pcEndRot: new Quaternion(),
}));

export function getSphereShape(r: number, margin?: number) {
    const sphere = new Ammo.btSphereShape(r);
    if (margin) sphere.setMargin(margin);
    return sphere;
}
    
export function getBoxShape(halfExtents: {x: number, y: number, z: number}, margin: number) {
    const ammoHalfExtents = helpers._ammoHalfExtents;
    ammoHalfExtents.setValue(halfExtents.x, halfExtents.y, halfExtents.z);
    const shape = new Ammo.btBoxShape(ammoHalfExtents);
    if (margin) shape.setMargin(margin);
    return shape;
}

export function getCapsuleShape(r: number, h: number, margin?: number) {
    const shape = new Ammo.btCapsuleShape(r, h);
    if (margin) shape.setMargin(margin);
    return shape;
}

export function getConvexHullShape(points: Readonly<Float32Array | number[]>) {
    // @ts-ignore
    const shape = new Ammo.btConvexHullShape(points, points.length, 3);
    return shape;
}

export function movePosition(rigidBody: RigidBodyComponent, newPosition: Vector3) {
    // TODO
    rigidBody.teleport(newPosition, rigidBody.entity.getRotation());
}

export function moveRotation(rigidBody: RigidBodyComponent, newRotation: Quaternion) {
    // TODO
    rigidBody.teleport(rigidBody.entity.getPosition(), newRotation);
}

export class RaycastHit extends RaycastResult {

    distance: number;

    constructor(entity: MaybeTransform, point: Vector3, normal: Vector3, hitFraction: number, distance: number) {
        super(entity, point, normal, hitFraction);
        this.distance = distance;
    }
    
    public clone() {

        return new RaycastHit(
            this.entity,
            Vec3Ex.clone(this.point),
            Vec3Ex.clone(this.normal),
            this.hitFraction,
            this.distance
        );
    }
}

export function getPerpendicularForward(up: Readonly<Vector3>, out = new Vector3()) {
    const forward = Vec3Ex.cross(DVector3.up, up, out);
    if (Vec3Ex.equals(forward, DVector3.zero)) {
        Vec3Ex.assign(forward, UnityVector3.forward);
    }
    return forward;
}

function assignSingle(pos?: Readonly<Vector3>, rot?: Readonly<Quaternion>) {

    if (pos) helpers._ammoPosFrom.setValue(pos.x, pos.y, pos.z);
    else     helpers._ammoPosFrom.setValue(0, 0, 0);

    if (rot) helpers._ammoRotFrom.setValue(rot.x, rot.y, rot.z, rot.w);
    else     helpers._ammoRotFrom.setValue(0, 0, 0, 1);

    helpers._ammoTransformFrom.setOrigin(helpers._ammoPosFrom);
    helpers._ammoTransformFrom.setRotation(helpers._ammoRotFrom);
}

function assignHelpersOrigin(from?: Readonly<Vector3>, to?: Readonly<Vector3>) {

    if (from) helpers._ammoPosFrom.setValue(from.x, from.y, from.z);
    else      helpers._ammoPosFrom.setValue(0, 0, 0);

    if (to) helpers._ammoPosTo.setValue(to.x, to.y, to.z);
    else    helpers._ammoPosTo.setValue(0, 0, 0);
}

function assignHelpersRotation(from?: Readonly<Quaternion>, to?: Readonly<Quaternion>) {
    
    if (from) helpers._ammoRotFrom.setValue(from.x, from.y, from.z, from.w);
    else      helpers._ammoRotFrom.setValue(0, 0, 0, 1);

    if (to) helpers._ammoRotTo.setValue(to.x, to.y, to.z, to.w);
    else    helpers._ammoRotTo.setValue(0, 0, 0, 1);
}

function assignHelpers(
    from?: Readonly<Vector3>,
    to?: Readonly<Vector3>,
    qFrom?: Readonly<Quaternion>,
    qTo?: Readonly<Quaternion>
) {
    assignHelpersOrigin(from, to);
    assignHelpersRotation(qFrom, qTo);

    helpers._ammoTransformFrom.setOrigin(helpers._ammoPosFrom);
    helpers._ammoTransformFrom.setRotation(helpers._ammoRotFrom);
    helpers._ammoTransformTo.setOrigin(helpers._ammoPosTo);
    helpers._ammoTransformTo.setRotation(helpers._ammoRotTo);
}

function convertHitsRayResults(
    convexCallback: Ammo.AllHitsRayResultCallback,
    results: RaycastHit[],
    maxDistance: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
) {
    let count = 0;

    if (convexCallback.hasHit()) {

        const collisionObjs = convexCallback.get_m_collisionObjects();
        const points = convexCallback.get_m_hitPointWorld();
        const normals = convexCallback.get_m_hitNormalWorld();
        const hitFractions = convexCallback.get_m_hitFractions();
        const numHits = Math.min(collisionObjs.size(), results.length);

        for (let i = 0; i < numHits; i++) {

            const obj  = collisionObjs.at(i);
            const body = Ammo.castObject(obj, Ammo.btRigidBody);

            if (body && body.entity) {

                if (filterTags && !body.entity.tags.has(...filterTags) ||
                    filterCallback && !filterCallback(body.entity)) {
                    continue;
                }

                const point = points.at(i);
                const normal = normals.at(i);
                const hitFraction = hitFractions.at(i);
                const distance = maxDistance * hitFraction;

                const raycastHit = 
                    //results[count] || (
                    results[count] = new RaycastHit(
                        body.entity,
                        new Vector3(point.x(), point.y(), point.z()),
                        new Vector3(normal.x(), normal.y(), normal.z()),
                        hitFraction,
                        distance
                    )//);

                raycastHit.entity = body.entity;
                raycastHit.hitFraction = hitFraction;
                raycastHit.distance = distance;
                //Vec3Ex.set(raycastHit.point, point.x(), point.y(), point.z());
                //Vec3Ex.set(raycastHit.normal, normal.x(), normal.y(), normal.z());

                //Debug.drawPoint({center: results[count].point});

                count++;
            }
        }
    }

    return count;
}

function convertHitsConvexResultsNonAlloc(
    convexCallback: Ammo.AllHitsConvexResultCallback,
    results: RaycastHit[],
    maxDistance: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
): number {

    let count = 0;

    if (convexCallback.hasHit()) {

        const collisionObjs = convexCallback.get_m_hitCollisionObjects();
        const points = convexCallback.get_m_hitPointWorld();
        const normals = convexCallback.get_m_hitNormalWorld();
        const hitFractions = convexCallback.get_m_hitFractions();
        const numHits = Math.min(collisionObjs.size(), results.length);

        for (let i = 0; i < numHits; i++) {

            const obj  = collisionObjs.at(i);
            const body = Ammo.castObject(obj, Ammo.btRigidBody);

            if (body && body.entity) {

                if (filterTags && !body.entity.tags.has(...filterTags) ||
                    filterCallback && !filterCallback(body.entity)) {
                    continue;
                }

                const point = points.at(i);
                const normal = normals.at(i);
                const hitFraction = hitFractions.at(i);
                const distance = maxDistance * hitFraction;

                const raycastHit =
                    //results[count] || (
                    results[count] = new RaycastHit(
                        body.entity,
                        new Vector3(point.x(), point.y(), point.z()),
                        new Vector3(normal.x(), normal.y(), normal.z()),
                        hitFraction,
                        distance
                    )//);
                
                raycastHit.entity = body.entity;
                raycastHit.hitFraction = hitFraction;
                raycastHit.distance = distance;
                //Vec3Ex.set(raycastHit.point, point.x(), point.y(), point.z());
                //Vec3Ex.set(raycastHit.normal, normal.x(), normal.y(), normal.z());

                //Debug.drawPoint({center: results[count].point});

                count++;
            }
        }
    }

    return count;
}

export function raycastFromToNonAlloc(
    system: PhysicsSystem,
    from: Readonly<Vector3>,
    to: Readonly<Vector3>,
    results: RaycastHit[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean
): number {

    //Debug.drawLine(from, to);

    assignHelpersOrigin(from, to);

    const maxDistance = Vec3Ex.distance(from, to);
    const rayCallback = new Ammo.AllHitsRayResultCallback(helpers._ammoPosFrom, helpers._ammoPosTo);
    const dynamicsWorld  = system.dynamicsWorld as Ammo.btDynamicsWorld;

    if (typeof filterCollisionMask === 'number') {
        rayCallback.set_m_collisionFilterMask(filterCollisionMask);
    }

    if (typeof filterCollisionGroup === 'number') {
        rayCallback.set_m_collisionFilterGroup(filterCollisionGroup);
    }

    dynamicsWorld.rayTest(helpers._ammoPosFrom, helpers._ammoPosTo, rayCallback);

    const count = convertHitsRayResults(rayCallback, results, maxDistance, filterTags, filterCallback);

    Ammo.destroy(rayCallback);

    return count;
}

export function raycastNonAlloc(
    system: PhysicsSystem,
    origin: Readonly<Vector3>,
    direction: Readonly<Vector3>,
    results: RaycastHit[],
    distance: number,
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean
) {
    const from  = origin;
    const to    = Vec3Ex.addRef2(from, Vec3Ex.multiplyScalarRef(Vec3Ex.normalize(direction, helpers._pcVec31), distance));
    
    return raycastFromToNonAlloc(system, from, to, results, filterCollisionMask, filterCollisionGroup, filterTags, filterCallback);
}

export function convexCastFromToNonAlloc(
    system: PhysicsSystem,
    fromPos: Readonly<Vector3>,
    toPos: Readonly<Vector3>,
    fromRot: Readonly<Quaternion>,
    toRot: Readonly<Quaternion>,
    shape: Ammo.btConvexShape,
    results: RaycastHit[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
) {
    assignHelpers(fromPos, toPos, fromRot, toRot);

    const maxDistance    = Vec3Ex.distance(fromPos, toPos);
    const convexCallback = new Ammo.AllHitsConvexResultCallback(helpers._ammoPosFrom, helpers._ammoPosTo);
    const dynamicsWorld  = system.dynamicsWorld as Ammo.btDynamicsWorld;

    if (typeof filterCollisionMask === 'number') {
        convexCallback.set_m_collisionFilterMask(filterCollisionMask);
    }

    if (typeof filterCollisionGroup === 'number') {
        convexCallback.set_m_collisionFilterGroup(filterCollisionGroup);
    }

    dynamicsWorld.convexSweepTest(shape, helpers._ammoTransformFrom, helpers._ammoTransformTo, convexCallback, allowedCcdPenetration || 0);

    const count = convertHitsConvexResultsNonAlloc(convexCallback, results, maxDistance, filterTags, filterCallback);
    
    Ammo.destroy(convexCallback);

    return count;
}

export function capsuleCastNonAlloc(
    system: PhysicsSystem,
    point1: Readonly<Vector3>,
    point2: Readonly<Vector3>,
    radius: number,
    direction: Vector3,
    results: RaycastHit[],
    maxDistance: number,
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
    bufShape?: Ammo.btCapsuleShape | Ammo.btCapsuleShapeX | Ammo.btCapsuleShapeZ
): number {

    const up       = Vec3Ex.subtract(point2, point1, helpers._pcVec31);
    const forward  = getPerpendicularForward(up, helpers._pcVec32);
    const rotation = QuatEx.lookRotation(forward, up, helpers._pcQuat1);
    const from     = Vec3Ex.lerp(point1, point2, 0.5, helpers._pcVec31);
    const to       = Vec3Ex.addRef2(from, Vec3Ex.multiplyScalarRef(Vec3Ex.normalize(direction, helpers._pcVec32), maxDistance));

    let shape = bufShape;

    if (!shape) {
        const height = Vec3Ex.distance(point1, point2);
        shape = getCapsuleShape(radius, height);
    }
    
    const count = convexCastFromToNonAlloc(
        system,
        from,
        to,
        rotation,
        rotation,
        shape,
        results,
        filterCollisionMask,
        filterCollisionGroup,
        filterTags,
        filterCallback,
        allowedCcdPenetration
    );

    if (!bufShape) {
        Ammo.destroy(shape);
    }

    return count;
}

export function sphereCastFromToNonAlloc(
    system: PhysicsSystem,
    from: Readonly<Vector3>,
    to: Readonly<Vector3>,
    radius: number,
    results: RaycastHit[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
    bufShape?: Ammo.btSphereShape
) {
    const capsuleShape = bufShape
        ? bufShape
        : getSphereShape(radius);
    
    const count = convexCastFromToNonAlloc(
        system,
        from,
        to,
        Quaternion.IDENTITY,
        Quaternion.IDENTITY,
        capsuleShape,
        results,
        filterCollisionMask,
        filterCollisionGroup,
        filterTags,
        filterCallback,
        allowedCcdPenetration
    );

    if (!bufShape) {
        Ammo.destroy(capsuleShape);
    }
    
    return count;
}

export function sphereCastNonAlloc(
    system: PhysicsSystem,
    origin: Readonly<Vector3>,
    radius: number,
    direction: Readonly<Vector3>,
    results: RaycastHit[],
    maxDistance: number,
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
    bufShape?: Ammo.btSphereShape
) {
    const from = origin;
    const to   = Vec3Ex.addRef2(from, Vec3Ex.multiplyScalarRef(Vec3Ex.normalize(direction, helpers._pcVec31), maxDistance));
    return sphereCastFromToNonAlloc(system, from, to, radius, results, filterCollisionMask, filterCollisionGroup, filterTags, filterCallback, allowedCcdPenetration, bufShape);
}

export function toShapePosition(collision: CollisionComponent, position: Readonly<Vector3>, rotation: Readonly<Quaternion>, hV: Vector3, hQ: Quaternion): Readonly<Vector3> {

    if (collision['_hasOffset']) {
        const lo  = collision.data.linearOffset;
        hQ.copy(rotation).transformVector(lo, hV);
        return hV.add(position);
    }

    return position;
}

export function toShapeRotation(collision: CollisionComponent, rotation: Readonly<Quaternion>, hQ: Quaternion): Readonly<Quaternion> {

    if (collision['_hasOffset']) {
        return hQ.copy(rotation).mul(collision.data.angularOffset);
    }

    return rotation;
}

export function computePenetration(
    system: PhysicsSystem,
    colliderA: pc.RigidBodyComponent,
    positionA: Readonly<Vector3>,
    rotationA: Readonly<Quaternion>,
    colliderB: pc.RigidBodyComponent,
    positionB: Readonly<Vector3>,
    rotationB: Readonly<Quaternion>,
    direction: OutObject<Vector3>,
    distance: Out<number | null>
) {
    const collisionA = colliderA.entity.collision!;
    const collisionB = colliderB.entity.collision!;
    const collisionAType = collisionA.type;
    const collisionAIsSimplest = collisionAType === 'capsule' || collisionAType === 'sphere' || collisionAType === 'box';

    // Simplest not support
    if (!collisionAIsSimplest) {
        Vec3Ex.assign(direction, DVector3.zero);
        distance.value = null;
        return false;
    }

    const dynamicsWorld = system.dynamicsWorld as Ammo.btDynamicsWorld;
    const colliderABody = colliderA.body as Ammo.btRigidBody;
    const colliderBBody = colliderB.body as Ammo.btRigidBody;

    const normalizePositionA = toShapePosition(collisionA, positionA, rotationA, helpers._pcVec31, helpers._pcQuat1);
    const normalizeRotationA = toShapeRotation(collisionA, rotationA, helpers._pcQuat1);
    const normalizePositionB = toShapePosition(collisionB, positionB, rotationB, helpers._pcVec32, helpers._pcQuat2);
    const normalizeRotationB = toShapeRotation(collisionB, rotationB, helpers._pcQuat2);

    assignHelpers(normalizePositionA, normalizePositionB, normalizeRotationA, normalizeRotationB);

    // Handle Gost A
    {
        helpers._ammoGhostObjectA.setCollisionShape(colliderABody.getCollisionShape());
        helpers._ammoGhostObjectA.setWorldTransform(helpers._ammoTransformFrom);

        if (helpers._ammoGhostObjectAWorld !== dynamicsWorld) {
            helpers._ammoGhostObjectAWorld = dynamicsWorld;
            helpers._ammoGhostObjectA.setCollisionFlags(pc.BODYFLAG_NORESPONSE_OBJECT);
            dynamicsWorld.addCollisionObject(helpers._ammoGhostObjectA, pc.BODYGROUP_NONE, pc.BODYMASK_NONE);
        }
    }

    // Handle Gost B
    {
        helpers._ammoGhostObjectB.setCollisionShape(colliderBBody.getCollisionShape());
        helpers._ammoGhostObjectB.setWorldTransform(helpers._ammoTransformTo);

        if (helpers._ammoGhostObjectBWorld !== dynamicsWorld) {
            helpers._ammoGhostObjectBWorld = dynamicsWorld;
            helpers._ammoGhostObjectB.setCollisionFlags(pc.BODYFLAG_NORESPONSE_OBJECT);
            dynamicsWorld.addCollisionObject(helpers._ammoGhostObjectB, pc.BODYGROUP_NONE, pc.BODYMASK_NONE);
        }
    }

    const contactResponse = new Ammo.ConcreteContactResultCallback();

    let minDepenetrationDistance  = 0;
    const minDepenetrationDirection = helpers._pcVec33;

    contactResponse.addSingleResult = (cp, colObj0Wrap) => {

        const manifold = Ammo.wrapPointer(cp, Ammo.btManifoldPoint);
        const wrapObj0 = Ammo.wrapPointer(colObj0Wrap, Ammo.btCollisionObjectWrapper);
        const colObj0  = wrapObj0.getCollisionObject();
        const ghostAIs0 = Ammo.compare(colObj0, helpers._ammoGhostObjectA);
        const directionSign = ghostAIs0 ? 1 : -1;
        const depenetrationDistance = manifold.getDistance();

        if (depenetrationDistance < 0 && depenetrationDistance < minDepenetrationDistance) {

            minDepenetrationDistance = depenetrationDistance;

            const ammoDirection = manifold.get_m_normalWorldOnB();

            Vec3Ex.set(minDepenetrationDirection, directionSign * ammoDirection.x(), directionSign * ammoDirection.y(), directionSign * ammoDirection.z());

            /*
            const ammoPointA    = manifold.get_m_positionWorldOnA();
            const ammoPointB    = manifold.get_m_positionWorldOnB();
            const pointA        = new Vector3(ammoPointA.x(), ammoPointA.y(), ammoPointA.z());
            const pointB        = new Vector3(ammoPointB.x(), ammoPointB.y(), ammoPointB.z());

            Debug.drawPoint({ center: pointA, color: pc.Color.GRAY });
            Debug.drawPoint({ center: pointB, color: pc.Color.CYAN });
            Debug.drawDirectionVector(minDepenetrationDirection, pointB);
            console.log(depenetrationDistance);
            */
            
            //const resultB = Vec3Ex.add(pointB, Vec3Ex.multiplyScalar(direction, collisionA.radius));

            //Debug.drawPoint({ center: resultB, color: pc.Color.BLACK });
            //Debug.drawPoint({ center: resultEntityB, color: pc.Color.MAGENTA });
        }

        return 1;
    }

    dynamicsWorld.contactPairTest(helpers._ammoGhostObjectA, helpers._ammoGhostObjectB, contactResponse);
    helpers._ammoGhostObjectA.setCollisionShape(helpers._ammoBlankSphere);
    helpers._ammoGhostObjectB.setCollisionShape(helpers._ammoBlankSphere);

    Ammo.destroy(contactResponse);

    if (minDepenetrationDistance !== 0) {

        //console.log(colliderB.entity.name, minDepenetrationDistance, minDepenetrationDirection.toString());

        Vec3Ex.assign(direction, minDepenetrationDirection);
        distance.value = Math.abs(minDepenetrationDistance);

        //console.log([colliderA.entity.getGuid(), colliderB.entity.getGuid()]);

        return true;
    }

    Vec3Ex.assign(direction, DVector3.zero);
    distance.value = null;
    return false;
}

const overlapEntities = new Set<MaybeTransform>();

export function overlapShapeNonAlloc(
    system: PhysicsSystem,
    shape: Ammo.btCollisionShape,
    position: Readonly<Vector3>,
    rotation: Readonly<Quaternion>,
    results: MaybeTransform[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number
) {
    const max = results.length - 1;
    const dynamicsWorld = system.dynamicsWorld as Ammo.btDynamicsWorld;

    assignSingle(position, rotation);

    helpers._ammoGhostObjectA.setCollisionShape(shape);
    helpers._ammoGhostObjectA.setWorldTransform(helpers._ammoTransformFrom);
    helpers._ammoOverlapContactResponse.addSingleResult = (
        cp: Ammo.btManifoldPoint,
        colObj0Wrap: Ammo.btCollisionObjectWrapper,
        partId0: number,
        index0: number,
        colObj1Wrap: Ammo.btCollisionObjectWrapper,
        partId1: number,
        index1: number
    ) => {

        if (overlapEntities.size > max) {
            return 0;
        }

        const manifold  = Ammo.wrapPointer(cp, Ammo.btManifoldPoint);
        const distance  = manifold.getDistance();

        if (distance > 0) {
            return 1;
        }

        let wrapper = Ammo.wrapPointer(colObj1Wrap, Ammo.btCollisionObjectWrapper);
        let colObj  = wrapper.getCollisionObject();

        const ghostIs1  = Ammo.compare(colObj, helpers._ammoGhostObjectA);
        if   (ghostIs1) {
            wrapper = Ammo.wrapPointer(colObj0Wrap, Ammo.btCollisionObjectWrapper);
            colObj  = wrapper.getCollisionObject();
        }

        const rigidBody = Ammo.castObject(colObj, Ammo.btRigidBody);

        if (!rigidBody || !rigidBody.entity) {
            return 1;
        }
    
        if (filterTags && !rigidBody.entity.tags.has(...filterTags) ||
            filterCallback && !filterCallback(rigidBody.entity)) {
            return 1;
        }
    
        overlapEntities.add(rigidBody.entity);

        /*
            const manifold = Ammo.wrapPointer(cp, Ammo.btManifoldPoint);
            const ammoDirection = manifold.get_m_normalWorldOnB();            
            const ammoPointA    = manifold.get_m_positionWorldOnA();
            const ammoPointB    = manifold.get_m_positionWorldOnB();
            const pointA        = new Vector3(ammoPointA.x(), ammoPointA.y(), ammoPointA.z());
            const pointB        = new Vector3(ammoPointB.x(), ammoPointB.y(), ammoPointB.z());
        
            Debug.drawPoint({ center: pointA, color: pc.Color.GRAY });
            Debug.drawPoint({ center: pointB, color: pc.Color.CYAN });
            Debug.drawDirectionVector(new Vector3(ammoDirection.x(), ammoDirection.y(), ammoDirection.z()), pointB);
        
            //console.log(body.entity.name);
        //*/

        return 1;
    }

    dynamicsWorld.contactTest(helpers._ammoGhostObjectA, helpers._ammoOverlapContactResponse);
    helpers._ammoGhostObjectA.setCollisionShape(helpers._ammoBlankSphere);

    let i = 0;
    for (let item of overlapEntities) results[i++] = item;

    overlapEntities.clear();

    return i;
}

export function overlapCapsuleNonAlloc(
    system: PhysicsSystem,
    point0: Readonly<Vector3>,
    point1: Readonly<Vector3>,
    radius: number,
    results: MaybeTransform[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
    bufShape?: Ammo.btCapsuleShape | Ammo.btCapsuleShapeX | Ammo.btCapsuleShapeZ,
): number {

    const capsuleShape = bufShape ? bufShape : getCapsuleShape(radius, Vec3Ex.distance(point1, point0));
    const up       = Vec3Ex.subtract(point1, point0, helpers._pcVec31);
    const forward  = getPerpendicularForward(up, helpers._pcVec32);
    const center   = Vec3Ex.lerp(point0, point1, 0.5, helpers._pcVec33);
    const rotation = QuatEx.lookRotation(forward, up, helpers._pcQuat1);

    const count = overlapShapeNonAlloc(
        system,
        capsuleShape, center, rotation,
        results,
        filterCollisionMask, filterCollisionGroup,
        filterTags, filterCallback,
        allowedCcdPenetration
    );

    if (!bufShape) {
        Ammo.destroy(capsuleShape);
    }

    return count;
}

export interface IPenetrationData {
    entity: MaybeTransform,
    normal: Vector3,
    distance: number,
}

const overlapAndComputePenetrationEntities = new Map<MaybeTransform, IPenetrationData>();

export function overlapShapeAndComputePenetrationNonAlloc(
    system: PhysicsSystem,
    shape: Ammo.btCollisionShape,
    position: Readonly<Vector3>,
    rotation: Readonly<Quaternion>,
    results: IPenetrationData[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
) {
    assignSingle(position, rotation);

    helpers._ammoGhostObjectA.setCollisionShape(shape);
    helpers._ammoGhostObjectA.setWorldTransform(helpers._ammoTransformFrom);
    helpers._ammoOverlapContactResponse.addSingleResult = (
        cp: Ammo.btManifoldPoint,
        colObj0Wrap: Ammo.btCollisionObjectWrapper,
        partId0: number,
        index0: number,
        colObj1Wrap: Ammo.btCollisionObjectWrapper,
        partId1: number,
        index1: number
    ) => {

        let wrapper  = Ammo.wrapPointer(colObj1Wrap, Ammo.btCollisionObjectWrapper);
        let colObj   = wrapper.getCollisionObject();
        let directionSign = 1;

        const ghostIs1 = Ammo.compare(colObj, helpers._ammoGhostObjectA);
        if (ghostIs1) {
            wrapper  = Ammo.wrapPointer(colObj0Wrap, Ammo.btCollisionObjectWrapper);
            colObj   = wrapper.getCollisionObject();
            directionSign = -1;
        }

        const rigidBody = Ammo.castObject(colObj, Ammo.btRigidBody);

        if (!rigidBody || !rigidBody.entity) {
            return 1;
        }
    
        if ((filterTags && !rigidBody.entity.tags.has(...filterTags)) ||
            (filterCallback && !filterCallback(rigidBody.entity))) {
            return 1;
        }

        // Ammo return negative normal for compound rigidbody
        const manifold = Ammo.wrapPointer(cp, Ammo.btManifoldPoint);
        const distance = manifold.getDistance();
        //console.log(ghostIsA, rigidBody.entity.name, distance, `[${-ammoDirection.x()}, ${-ammoDirection.y()}, ${-ammoDirection.z()}]`);

        if (distance < 0) {

            const normalizeDistance = Math.abs(distance);

            let data = overlapAndComputePenetrationEntities.get(rigidBody.entity);
            if (data) {
                if (data.distance > normalizeDistance) {
                    data.distance = normalizeDistance;

                    const ammoDirection = manifold.get_m_normalWorldOnB();

                    Vec3Ex.set(data.normal, directionSign * ammoDirection.x(),  directionSign * ammoDirection.y(),  directionSign * ammoDirection.z());
                }
            }
            else {
                const ammoDirection = manifold.get_m_normalWorldOnB();
                data = {
                    entity:   rigidBody.entity,
                    distance: normalizeDistance,
                    normal:   new Vector3(directionSign * ammoDirection.x(),  directionSign * ammoDirection.y(),  directionSign * ammoDirection.z())
                }
                overlapAndComputePenetrationEntities.set(rigidBody.entity, data);
            }
        }

        /*
        const ammoDirection = manifold.get_m_normalWorldOnB();
        const ammoPointA = manifold.get_m_positionWorldOnA();
        const ammoPointB = manifold.get_m_positionWorldOnB();
        const pointA     = new Vector3(ammoPointA.x(), ammoPointA.y(), ammoPointA.z());
        const pointB     = new Vector3(ammoPointB.x(), ammoPointB.y(), ammoPointB.z());

        Debug.drawPoint({ center: pointA, color: pc.Color.GRAY });
        Debug.drawPoint({ center: pointB, color: pc.Color.CYAN });
        Debug.drawDirectionVector(new Vector3(ammoDirection.x(), ammoDirection.y(), ammoDirection.z()), pointB);
        */
        
        return 1;
    }

    const dynamicsWorld = system.dynamicsWorld as Ammo.btDynamicsWorld;

    dynamicsWorld.contactTest(helpers._ammoGhostObjectA, helpers._ammoOverlapContactResponse);
    helpers._ammoGhostObjectA.setCollisionShape(helpers._ammoBlankSphere);

    let i = 0;
    for (let item of overlapAndComputePenetrationEntities) {
        results[i++] = item[1];
    }

    overlapAndComputePenetrationEntities.clear();

    return i;
}

export function overlapCapsuleAndComputePenetrationNonAlloc(
    system: PhysicsSystem,
    point0: Readonly<Vector3>,
    point1: Readonly<Vector3>,
    radius: number,
    results: IPenetrationData[],
    filterCollisionMask?: number,
    filterCollisionGroup?: number,
    filterTags?: Readonly<string[]>,
    filterCallback?: (entity: pc.Entity) => boolean,
    allowedCcdPenetration?: number,
    bufShape?: Ammo.btCapsuleShape | Ammo.btCapsuleShapeX | Ammo.btCapsuleShapeZ,
) {
    const capsuleShape = bufShape ? bufShape : getCapsuleShape(radius, Vec3Ex.distance(point1, point0));
    const up       = Vec3Ex.subtract(point1, point0, helpers._pcVec31);
    const forward  = getPerpendicularForward(up, helpers._pcVec32);
    const center   = Vec3Ex.lerp(point0, point1, 0.5, helpers._pcVec33);
    const rotation = QuatEx.lookRotation(forward, up, helpers._pcQuat1);

    const count = overlapShapeAndComputePenetrationNonAlloc(
        system,
        capsuleShape,
        center,
        rotation,
        results,
        filterCollisionMask,
        filterCollisionGroup,
        filterTags,
        filterCallback,
        allowedCcdPenetration
    );

    if (!bufShape) {
        Ammo.destroy(capsuleShape);
    }

    return count;
}

export const Physics = {
    movePosition,
    moveRotation,

    overlapShapeNonAlloc,
    overlapShapeAndComputePenetrationNonAlloc,

    // Unity friendly
    convexCastFromToNonAlloc,
    raycastNonAlloc,
    raycastFromToNonAlloc,
    overlapCapsuleNonAlloc,
    capsuleCastNonAlloc,
    sphereCastNonAlloc,
    sphereCastFromToNonAlloc,
    computePenetration,
    overlapCapsuleAndComputePenetrationNonAlloc,
}

export default Physics;
