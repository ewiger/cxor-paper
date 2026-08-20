# cxor

Context Transformation Algorithms (CTA) denotes a broader family of context-addressed transforms.
This repository studies one deliberately narrow binary member, **cxor** (CTA-XOR): contexts of
order $m$, one bit of state per context, an XOR residual emitter, and a one-bit update rule.

For an order-$m$ binary model, the previous $m$ source bits select a one-bit table entry
$\hat{x}$. cxor emits

$$
r=x\oplus\hat{x},
$$

so `0` records a correct prediction and `1` an error. The earlier implementation used
$\mathsf{W}_5(s,x)=x$: the next time a context occurs, predict the successor observed on its
previous occurrence. cxor preserves length and does not compress by itself; its possible role is
to expose contextual predictability for a downstream coder.

Allowing the one-bit update to be any Boolean function of the old cell and current symbol gives
exactly sixteen wirings within this cxor model. The paper classifies their causal invertibility,
per-context operators, and adaptation horizons, with $\mathsf{W}_5$ as the historical protagonist and
$\mathsf{W}_6(s,x)=s\oplus x$ as its running-parity counterpart.

```bash
make            # build cxor-paper.pdf
make check      # run the verification harness
```

The historical implementation is [`cta.py`](cta.py), and the complete one-bit verification harness is
[`cta_spec.py`](cta_spec.py).

Author: Yauhen Yakimovich
