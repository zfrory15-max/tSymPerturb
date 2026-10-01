# tSymPerturb manuscript method contract

Source: Zheng Zhu et al., *tSymPerturb converts longitudinal symptom networks into time-indexed intervention strategies*, user-supplied 22-page manuscript, arXiv submission header dated 17 August 2026. Page references below are the printed PDF pages. The attachment has a `.doc` suffix but is a PDF. No supplementary appendix is present in this source: Methods are section 9, pp. 12–20; references occupy pp. 21–22. Equations 13–53 were cross-checked against rendered pages 13–18. This is a source-grounded contract, not a claim to reproduce the original unpublished code or synthetic population.

## 1. Objects, scales and interpretation

All p symptoms have the same ordered identities at source and outcome and are oriented so higher is worse. B is p by p, with B[j,i] the path FROM source i TO outcome j. Columns, never rows, are outgoing profiles. Diagonal entries are autoregressive paths. Responses are baseline minus perturbed predictions; positive is improvement, negative is predicted worsening. Never silently discard negative responses.

The CLPN is fitted once for a given dataset. Query the fitted model for hypothetical perturbations; do not refit after each perturbation. Re-estimation occurs across independent datasets or participant bootstrap replicates, not inside each virtual intervention. An exact constant target causes no singularity in a forward query. No matrix inverse or cross-sectional equilibrium operation belongs in this reference method.

Fitted-scale B, source means, anchors, covariance, residual covariance and outcome standard deviations must share a coherent coordinate system. If standardizing input and outcome separately, transform anchors on the original source scale first: c_z=(c_raw-mu_source_raw)/sd_source_raw (Eq. 12). Standardized zero represents a sample mean, not symptom absence. Use the UNPERTURBED follow-up SD as the fixed denominator in standardized response (Eqs. 22,46), not the perturbed SD.

## 2. Equation inventory and implementation meanings

### Model and source moments (pp. 3,12–14)

- Eqs. 1 and 8: X2=a+B X1+epsilon; E(epsilon)=0; Var(epsilon)=Psi. Eq. 8 adds Cov(X1,epsilon)=0.
- Eq. 9: mu2=a+B mu1. Eq. 10: Sigma2=B Sigma1 B^T+Psi. Eq. 11: Cov(X1,X2)=Sigma1 B^T.
- Eq. 12: c_i,z=(c_i-mu1_i)/sigma1_i.
- Eqs. 2 and 13: for target set S, X1,S*(d)=cS+D_mu(d)(mu1,S-cS)+D_sigma(d)(X1,S-mu1,S). D maps are diagonal; non-target source variables are unchanged. D_mu(0)=D_sigma(0)=I; the exact anchor endpoint uses D_mu(1)=D_sigma(1)=0.
- Eq. 14: mu1,S*=cS+D_mu(mu1,S-cS).
- Eq. 15: Sigma1,SS*=D_sigma Sigma1,SS D_sigma^T.
- Eq. 16: in non-target,target order (K,S), covariance blocks are Sigma_KK, Sigma_KS D_sigma^T, D_sigma Sigma_SK, D_sigma Sigma_SS D_sigma^T. Equivalent implementation: full diagonal A has scale factors on S and ones on K; Sigma1*=A Sigma1 A^T. Scale cross-covariances too.
- Eqs. 17 and 18: mu2*=a+B mu1*; Sigma2*=B Sigma1* B^T+Psi. The intercept and residual covariance stay fixed.
- Eqs. 3 and 19: R_S=mu2-mu2*=B(mu1-mu1*).
- Eq. 20: exact t-vKO sets target mean to its anchor and target variance to zero, with associated target cross-covariances zero.
- Eqs. 4 and 21: R_i^vKO=B[:,i] delta_i, where scalar delta_i=mu1_i-c_i.
- Eq. 22: Delta[j<-i]=B[j,i] delta_i / s2_j; s2_j=sqrt(Sigma2[j,j]).
- Eq. 23: linked t-vKD X1,i*=c_i+(1-d)(X1,i-c_i), d in (0,1); support d=0 and 1 as baseline and t-vKO endpoints for testing/dosage curves.
- Eq. 24: target mean c_i+(1-d)(mu1_i-c_i); target variance (1-d)^2 sigma1_i^2.
- Eqs. 5 and 25: R_i(d)=d B[:,i] delta_i.
- Eq. 26: for a fixed linear weighted utility of standardized means, G_i(d)=d G_i(1).
- Eq. 27: G_i(d)/d=G_i(1) for d>0 and derivative at zero equals G_i(1). Do not divide by zero at d=0. Efficiency and low-dose responsiveness add no separate score dimensions in this reference case.
- Eq. 28: clipped Gaussian expectation for Y~N(m,s^2), h(y)=min(U,max(L,y)), a=(L-m)/s, b=(U-m)/s:
  E[h(Y)]=L Phi(a)+m[Phi(b)-Phi(a)]+s[phi(a)-phi(b)]+U[1-Phi(b)]. At s=0 use h(m). This is a clipped/winsorized normal, NOT the conditional expectation of a truncated normal. Covariance propagation is needed to obtain perturbed s. Clipping E[Y] alone is generally incorrect. This optional boundary model may produce curvature and nonadditivity; label it explicitly and do not apply the unbounded nulls as if clipping were absent.

### Transition operators (pp. 5,15)

- Eqs. 6 and 29: edge block B*=B-q B[j,i] e_j e_i^T, q in [0,1]. Only B[j,i] is multiplied by 1-q.
- Eq. 30: at prespecified reference source profile x_ref, signed predicted difference is (B-B*)x_ref=q B[j,i] x_ref_i e_j.
- Eq. 31: average source-distribution effect is q B[j,i] mu1_i e_j. If centered, this is zero; do not interpret coefficient magnitude as a signed benefit.
- Eq. 32: full source-node block B*=B-q B[:,i] e_i^T; attenuate the whole outgoing column including B[i,i].
- Eq. 33: cross-only variant B*=B-q(B[:,i]-B[i,i]e_i)e_i^T; preserve autoregressive diagonal. It is a supported alternative, not the full-block reference score.
- Eq. 34: Q_H(B)=1^T [sum over h=1..H of (gamma abs(B))^h] 1, 0<gamma<1. Absolute value is elementwise BEFORE matrix powers. Powers are matrix products, not elementwise powers; include h=1 through H and retain diagonals. Finite horizon avoids requiring an inverse or infinite-series convergence.
- Eq. 35: R_comm_i=[Q_H(B)-Q_H(B_block_i)]/Q_H(B). Table 2 (p.6) specifies q=.8 for the reference tVPPS. H and gamma are not numerically specified; declare configured choices. This unsigned quantity is route capacity, not signed clinical efficacy. Keep signed reference-state effects available alongside it.

### Three strategy procedures (pp. 5–6,14–17)

The FOUR primitives are t-vKO, t-vKD, edge communication block and source-node communication block. t-vDP is a dose-response evaluation procedure; combination construction and sequence optimization are the other TWO procedures. Do not describe seven independent primitive operators.

1. **t-vDP** evaluates a prespecified dose grid with the chosen state map and outcome functional. Linear linked reference responses must satisfy Eq. 26; thresholds, boundaries, nonlinear maps, interactions or state-dependent B must be declared if introduced.
2. **Combination construction:** Eq. 36: R_{i,k}(d_i,d_k)=B(d_i delta_i e_i+d_k delta_k e_k)=R_i(d_i)+R_k(d_k). Eq. 37: on the SAME fixed outcome comparison set with identical weights/SDs, G_ik=G_i+G_k and NA_ik=G_ik-G_i-G_k=0. Independent source changes add even if source symptoms were correlated.
   - Eq. 38: I_ik=G_ik^(-ik)-max(G_i^(-ik),G_k^(-ik)). ALL three utilities exclude BOTH outcomes i and k. Do not compare each single's own cross-only set or compare full-outcome pairs with target-excluding singles.
   - Eq. 39: R_comb_i=(1/|T_i|) sum_{k in T_i} max(0,I_ik), with prespecified partner set T_i. Preserve signed I_ik in output. Positive incremental pair value is not synergy; for the linear case on the same set it equals min(G_i^(-ik),G_k^(-ik)). Partner sets must contain distinct valid non-self partners.
3. **Sequence:** distinguish the following two mathematically different problems:
   - Eq. 40: observed multi-wave X_{t+1}=a_t+B_t X_t+epsilon_t.
   - Eq. 41: for a one-time source mean improvement vector delta_t, R_{t+h}=B_{t+h-1}...B_t delta_t. Apply matrices in chronological order to the vector, never reverse them. At h=1 this is B_t delta_t.
   - Eq. 42: stationary extension R_{t+h}=B^h delta_t. Reusing one two-wave B is an assumption, not additional longitudinal evidence.
   - Eq. 43: repeated intervention recursion Delta^-_{t+1}=B_t(Delta^-_t+u_t), where Delta^-_t=mu_t^0-mu_t^- and u_t is the additional source improvement immediately BEFORE transition t. For a state-dependent anchoring rule, compute u_t from the currently perturbed state, not the baseline mean. The manuscript does not specify one universal dose/action rule or optimization algorithm.
   - Eq. 44: finite-horizon utility U_i^(H)=sum_{h=1..H} eta^(h-1) [w^T R_{t+h,i}/sum_j w_j], 0<eta<=1. As printed this uses R (raw-scale improvement), not Delta (SD-standardized). Name/report the chosen scale and do not silently change the equation. The manuscript notation R_{t+h,i} denotes the response profile for intervention i.
   - Eq. 45: TWO-WAVE target-acquisition order pi, prefix sets S_m, J(pi)=sum_{m=1..L} gamma^(m-1)[U(S_m)-U(S_{m-1})]-sum_{m=1..L} lambda_{pi_m}. This orders decisions to add targets, not biological treatment timing. For additive U, fixed final set, gamma=1 and order-independent costs, J is invariant to permutation. With gamma<1 the discounted objective can differ by acquisition order even for additive U; do not claim blanket invariance under discounting. This acquisition gamma is conceptually distinct from Q_H's propagation discount.

### Seven reference utilities and tVPPS (pp. 6,17–18)

All source-state utilities used in the reference score evaluate d=1. Generic efficacy functions may also evaluate partial doses.

- Eq. 46: Delta[j<-i](d)=[mu2_j-mu2_j^(i,d)]/s2_j.
- Eq. 47 downstream efficacy: weighted mean over ALL p outcomes, sum_j w_j Delta_j / sum_j w_j.
- Eq. 48 cross-symptom spillover: weighted mean over j != i, with denominator also restricted to j != i.
- Eq. 49 breadth: unweighted fraction of p-1 non-target outcomes with Delta_j >= tau; Table 2 reference tau=.05 SD. Inclusive comparison.
- Eq. 50 cross-module reach: among M-1 modules OTHER than target's module, fraction whose UNWEIGHTED within-module mean Delta >= tau_m; Table 2 reference tau_m=.03 SD. First average each module, then count modules. Do not pool symptoms or weight modules by size.
- Eq. 35 communication-block value as above, q=.8 and explicit H/gamma.
- Eq. 39 combination value as above, with explicit partner sets and common comparison sets.
- Eq. 51 positive spillover fraction: sum_{j!=i} w_j max(Delta_j,0) / [sum_j w_j max(Delta_j,0)+1e-12]. Keep manuscript epsilon. It is a fraction of POSITIVE response, not signed net benefit. All-zero/nonpositive improvement produces zero with nonnegative weights.
- Eq. 52: direction-align each raw utility, then candidate-set min-max score=100*(R_im-min_l R_lm)/(max_l R_lm-min_l R_lm). A constant dimension scores 50 for every candidate. All seven reference definitions are larger-is-higher utility by construction, with the stated limitations. Normalize over the declared target candidate set, not all symptom nodes unless all are candidates.
- Eqs. 7 and 53: tVPPS_i=sum_{m=1..7} omega_m score_im / sum_m omega_m, omega_m>=0. Require strictly positive total weight. Equal weights are a methodological verification choice, not validated clinical utilities. Report raw seven-dimensional profiles next to the composite. Uncertainty is NOT an eighth dimension; dose efficiency/low-dose responsiveness are NOT extra dimensions.

### Verification, estimation and uncertainty (pp. 7,9,18–20)

- Eq. 54: simulation intercept a=(I-B)mu1 to keep unperturbed follow-up mean equal to baseline mean. This is a generating-system choice, not required of an applied CLPN.
- Eq. 55: Monte Carlo response estimate is sample baseline X2 mean minus sample perturbed X2 mean, then divide by unperturbed s2_j for standardized comparisons. Draws verify expectations; they are not independent clinical datasets.
- Eq. 56: absolute outgoing-strength benchmark sum_j |B[j,i]| (includes diagonal).
- Eq. 57: signed outgoing expected-influence benchmark sum_j B[j,i] (includes diagonal).
- Eq. 58: top-K bootstrap probability P_i^(K)=mean_b indicator(rank_i^(b)<=K). Bootstrap the entire participant-level pipeline: paired wave resampling, standardization, anchor transformation, CLPN fitting, perturbation, utility calculation, candidate-set normalization and ranking. Merely perturbing a single fitted B is not complete-pipeline uncertainty.

The verification used 22 synthetic nodes in four modules, 53 nonzero off-diagonal paths (47 positive, 6 negative), autoregressive coefficients .43–.71, spectral radius .827, means 1.44–2.64, source SDs .63–.83, residual SDs .38–.52, and master seed 20260816. It fit each outcome with lasso alpha=.03 without intercept after per-wave standardization. There were 200 independently generated datasets each at n=250,500,1000, plus 250,000 Monte Carlo draws per target. These are facts about the manuscript verification, not a universal estimator, validated sample-size rule, or defaults required of all applications. Reported numerical results cannot be recreated from ranges and a seed without the actual generating parameters and algorithms.

## 3. Required invariants and implementation tests

1. Dose zero recovers baseline moments and zero response; dose one recovers the exact anchor and zero targeted source variance/covariance.
2. An asymmetric B fixture verifies columns are outgoing, and edge i->j changes exactly B[j,i].
3. Linked dose d gives R(d)=d R(1); fixed linear utility gives G(d)=d G(1). Test signed coefficients and unequal SDs.
4. Joint state response equals the sum of singles; NA=0 on a common outcome set. Pair exclusion and denominators are independently checked.
5. Full node block attenuates the diagonal; cross-only leaves it unchanged. q=0 recovers B; q=1 zeros the chosen coefficient/column subset.
6. Q_H is a finite sum of matrix powers of gamma*abs(B). Negative-edge fixtures distinguish abs-before-power from abs-after-power; off-diagonal fixtures distinguish matrix from elementwise powers.
7. Source covariance scaling affects target cross-covariances, not non-target covariance. Covariance remains symmetric PSD within numerical tolerance.
8. Unequal-size modules verify per-module unweighted mean and equal module counting. A response exactly at a threshold counts.
9. Constant score dimensions give 50; normalized nonconstant dimensions range 0–100; a singleton candidate set gives all-neutral dimensions, not false precision.
10. Noncommuting wave matrices verify chronological products; repeated intervention is applied before transition; no-intervention recursion matches the one-time response.
11. Clipped-normal tests cover s=0, no effective clipping, symmetric zero mean, and far-tail bounds. Reference linear nulls remain separate from bounded-mode checks.
12. Fixed-set additive acquisition with no discount is order invariant; discounted acquisition and repeated-wave intervention must not be conflated.
