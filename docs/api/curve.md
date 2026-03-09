# API Reference: qaws_curve.h

Base curve lifecycle.

---

## qaws_curve_destroy

```c
void qaws_curve_destroy(qaws_curve* curve);
```

Destroys a curve object and frees all associated memory.

**Parameters:**
- `curve` -- Curve to destroy. May be NULL (no-op).

**Semantics:**
- Must be called for every curve created with `qaws_curve_create_*`.
- After destruction, the pointer must not be used.
- Any traversal objects referencing this curve become invalid.
