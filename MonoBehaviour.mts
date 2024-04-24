import { MaybeTransform } from "./Types.mjs";

export interface IClass<T> {
    new (...args: any[]): T;
}

export type Class<T> = IClass<T> | (Function & { prototype: T});

export type PlayCanvasScript = pc.ScriptType;
export const PlayCanvasScript = pc.ScriptType;

export function getComponentByName<TScript extends PlayCanvasScript = PlayCanvasScript>(entity: MaybeTransform, name: string): TScript | null {
    if (!entity.script)
        return null;
    return entity.script.get(name) as unknown as TScript;
}

export function getComponentByType<TScript extends PlayCanvasScript = PlayCanvasScript>(entity: MaybeTransform, type: Class<TScript>): TScript | null {

    if (!entity.script)
        return null;
    
    return entity.script.get(type as any) as unknown as TScript;
}

export function getFirstParentsComponentByType<TScript extends PlayCanvasScript = PlayCanvasScript>(entity: MaybeTransform, type: Class<TScript>, ascentCount: number = 3): TScript | null {

    let ascentEntity = entity.parent;

    for (let i = 0; i < ascentCount; i++) {

        if (!ascentEntity) {
            return null;
        }
        
        const component = getComponentByType(ascentEntity as MaybeTransform, type);

        if (component) {
            return component;
        }

        ascentEntity = ascentEntity.parent;
    }
    
    return null;
}

export abstract class MonoBehaviour<T extends MonoBehaviour<T> = MonoBehaviour<any>> extends PlayCanvasScript<T> {

    private ____monoBehaviourDidStart: boolean = false;

    readonly initialize = (): void => {

        this.on('destroy', this.____monoBehaviourDestroy, this);
        this.on('state', this.____monoBehaviourOnStateChange, this);
        this.on('attr', this.____monoBehaviourOnChangeAttr, this);

        this.awake();

        if (this.enabled) {
            this.onEnable();
        }
    }

    readonly postInitialize = (): void => {
        if (this.enabled) {
            this.____monoBehaviourDidStart = true;
            this.start();
        }
    }

    private ____monoBehaviourDestroy() {
        this.off('destroy', this.____monoBehaviourDestroy, this);
        this.off('state', this.____monoBehaviourOnStateChange, this);
        this.off('attr', this.____monoBehaviourOnStateChange, this);
        this.onDestroy();
    }

    private ____monoBehaviourOnChangeAttr<TKey extends keyof T = any>(name: TKey, value: any, prevValue: any) {

        return this.onChangeAttribute(name, value, prevValue);
    }

    private ____monoBehaviourOnStateChange(state: boolean) {

        if (state && !this.____monoBehaviourDidStart) {
            this.____monoBehaviourDidStart = true;
            this.start();
        }

        if(state) this.onEnable();
        else      this.onDisable();
    }

    /*
    readonly update = (dt: number): void => {
        //this.fixedUpdate(dt);
        this.update(dt);
    }
    */
    
    readonly postUpdate = (dt: number): void => {
        this.lateUpdate(dt);
    }

    public awake() {}
    public start() {}
    public onEnable() {}
    public onDisable() {}
    public onDestroy() {}

    //public fixedUpdate(dt: Array<string>) {}
    //public update(dt: number) {}
    public lateUpdate(dt: number) {}
    
    public onChangeAttribute<TKey extends keyof T>(attribute: TKey, value: any, prev: any): void {}

    public getComponent<TScript extends PlayCanvasScript = PlayCanvasScript>(name: string): TScript | null;
    public getComponent<TScript extends PlayCanvasScript>(type: Class<TScript>): TScript | null;
    public getComponent(nameOrType: any): any {
        return this.entity.script ? this.entity.script.get(nameOrType) : null;
    }
    
    public getComponentOrThrow<TScript extends PlayCanvasScript = PlayCanvasScript, TError = undefined>(name: string, error?: TError): TScript;
    public getComponentOrThrow<TScript extends PlayCanvasScript, TError = undefined>(type: Class<TScript>, error?: TError): TScript;
    public getComponentOrThrow(nameOrType: any, error?: any): any {
        
        let script: any
        if (!this.entity.script || !(script = this.entity.script.get(nameOrType))) {
            throw error || new Error('Script not exists');
        }

        return script;
    }
}

export default MonoBehaviour;
