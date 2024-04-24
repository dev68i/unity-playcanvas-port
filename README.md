This functionality makes it easier to transfer scripts from Unity to Play Canvas.
Let me give you a few examples:

```
using UnityEngine;

public class RotatePlayer : MonoBehaviour
{
    void Update()
    {
        if (Input.GetAxis("Horizontal") > 0)
        {
            Vector3 rotate = transform.eulerAngles;
            rotate.y = 0;
            transform.rotation = Quaternion.Euler(rotate);
        }
        
        if (Input.GetAxis("Horizontal") < 0)
        {
            Vector3 rotate = transform.eulerAngles;
            rotate.y = 180;
            transform.rotation = Quaternion.Euler(rotate);
        }
    }
}
```

Transfer to:
```
import MonoBehaviour from './MonoBehaviour.mjs';
import QuternionEx from './QuaternionEx.mjs';

export class RotatePlayer extends MonoBehaviour {

    update() {

        if (this.app.keyboard.isPressed(pc.KEY_A)) {
            const rotate = QuternionEx.toEulerAngles(this.entity.getRotation());
            rotate.y = 0;
            this.entity.setRotation(QuternionEx.euler(rotate));
        }

        if (this.app.keyboard.isPressed(pc.KEY_D)) {
            const rotate = QuternionEx.toEulerAngles(this.entity.getRotation());
            rotate.y = 180;
            this.entity.setRotation(QuternionEx.euler(rotate));
        }
    }
}
```
