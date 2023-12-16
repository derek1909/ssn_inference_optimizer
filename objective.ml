open Printf
open Math.Core
open Autodiff
open Util
module F = Operators.On_floats
module PE = Pervasives

(*
let () = Printexc.record_backtrace true
*)
(* Objective function defined in Autodiff, for learning a collection of moments.

   Free parameters:
   - 4 synaptic weights profiles ({E/I}-{E/I}), assumed to be symmetric and translation invariant
   - Cholesky factor (UPPER triangular) of the input noise covariance matrix
   - mu for the auxiliary neurons
   - Cholesky factor (UPPER triangular) involving the auxiliary neurons

*)

(* set of functions to find and report NANs or undesired neg*)

let is_nan x = 
  match classify_float x with FP_nan -> true | _ -> false

let vec_contains_nan v = 
  Vec.fold (fun accu x -> accu || is_nan x) false v

let report_nan x x_name =
  if (vec_contains_nan x) then printf "%s contains a NaN!\n%!" x_name

let is_zero x = 
  if (x = 0.0) then true else false  

let vec_contains_zero v = 
  Vec.fold (fun accu x -> accu || is_zero x) false v      

let report_zero x x_name =
  if (vec_contains_zero x) then printf "%s contains a 0!\n%!" x_name

let is_neg x = 
  if (x < 0.0) then true else false  

let vec_contains_neg v = 
  Vec.fold (fun accu x -> accu || is_neg x) false v      

let report_neg x x_name =
  if (vec_contains_neg x) then printf "%s contains a negative value!\n%!" x_name


(* ell is the UPPER triangular Cholesky factor of sigma *)
type target = { 
  h_vec: vec; (* vector of size M (number of sites on the ring) *)
  mu_vec: vec; (* M-vector *)
  sigma_mat: mat; (* MxM matrix *)
}

module Make (P: sig 
    (* Basic parameters *)
    val verbose: bool
    val m: int (* number of sites on the ring *)
    val k: float (* scaling constant of the powerlaw *)
    val tau_e: vec (* E cell membrane time constants *)
    val tau_i: vec (* I cell membrane time constants *)
    val tau_eta: float (* noise autocorrelation time *)
    val dt: float (* dt for the forward integration *)
    val t_max: float (* maximum integration time *)
    val t_subsamp: float (* time-span between points to integrate in the cost *) 
    val t_max_slow: float (* maximum integration time for the slowness cost *)
    val input_noise_size : float (* desired size of the diag. elements of sigma_eta *)
    val input_baseline_lb : float (* absolute minimum for the input baseline - otherwise you'll get nan *)
    val lambda_mean: float (* Lagrange multipliers *)
    val lambda_var: float
    val lambda_cov: float
    val lambda_slow: float
    val lambda_noise: float
    val lambda_noise_width: float
    val lambda_rates: float
    val lambda_sigma: float
    val lambda_reg_ratio: float
    val lambda_reg_var: float
    (* Targets *)
    val targets: target array
    (* sample-based estimate *)
    val n_trials: int option
  end) = struct

  open P
  let () = assert (m mod 2 = 0)
  let n = 2*m
  let half_m = m/2
  let half_m_plus_one = half_m + 1
  let n_targ = Array.length targets
  let taus = Vec.append tau_e tau_i
  let _ = assert (Vec.dim taus = n)
  let n_time_bins = round (t_max /. dt)
  let n_time_bins_slow = round (t_max_slow /. dt)
  let subsamp_bins = round (t_subsamp /. dt)
  let () = if verbose then printf "n_time_bins = %i\n%!" n_time_bins
  let min_time_bins = round (0.05 /. dt)
  let n_prms_base = 3 (* input_baseline, input_scaling, input_nl_pow *)
  let n_prms_w = 8 (* for each of 4 quadrants: height and width *)
  let n_prms_sigma_eta = 3 (* 2x2 covariance matrix, Kronecker-produced with... *) 
                         + 1 (* ... a base posdef matrix with squared exponential kernel *)
  let n_prms = n_prms_base + n_prms_w + n_prms_sigma_eta

  let () = if verbose then printf "Optimizing over %i free parameters\n%!" n_prms

  (* give the multipliers their natural scaling with system size *)
  let lambda_mean = lambda_mean  /. float (2*m*n_targ*n_time_bins/subsamp_bins)
  let lambda_var = lambda_var  /. float (2*m*n_targ*n_time_bins/subsamp_bins)
  let lambda_cov = lambda_cov  /. float (2*m*m*n_targ*n_time_bins/subsamp_bins)
  let lambda_slow = lambda_slow /. float (2*m*n_time_bins_slow*n_targ)
  let lambda_noise = lambda_noise /. float n
  let lambda_noise_width = lambda_noise_width
  let lambda_rates = lambda_rates /. float n
  let lambda_sigma = lambda_sigma /. float n
  let lambda_reg_ratio = lambda_reg_ratio /. float (m*n_targ)
  let lambda_reg_var = lambda_reg_var/. float (m*n_targ)
   /. float m
  (* Some useful constants *)

  let two_k = 2. *. k
  let normal_pdf_const = 1. /. sqrt (2. *. Gsl.Math.pi)
  let minus_inv_tau_eta = -. 1. /. tau_eta
  let inv_taus = Vec.map (fun tau -> 1. /. tau) taus
  let dt_inv_taus = Vec.map (fun tau -> dt /. tau) taus
  let eps1 = 1. -. dt /. tau_eta
  let eps2 = dt *. (1. +. dt /. tau_eta)
  let eps3 = sqrt (2.0 *. dt /. tau_eta)

  (* -----------------------------------------------------------------------
     @@   Parameter packing/unpacking                                     @@
     ----------------------------------------------------------------------- *)

  type 'a t_weight_prms = { width: 'a; height: 'a }

  type 'a t_weight_profile = { ee: 'a t_weight_prms;
                               ei: 'a t_weight_prms;
                               ie: 'a t_weight_prms;
                               ii: 'a t_weight_prms }

  type 'a t_sigma_eta_prms = { width: 'a; std_e: 'a; std_i: 'a; rho: 'a }

  type ('a, 'b, 'c) parameters = {
    input_baseline: 'a;
    input_scaling: 'a;
    input_nl_pow: 'a;
    weight_prms: 'a t_weight_profile;
    sigma_eta_prms: 'a t_sigma_eta_prms }

  (* unpack the big vector of parameters *)
  let%diff unpack p =
    let open O in
    let c = ref 0 in
    let pop () = incr c; V.get p !c in
    let input_baseline = input_baseline_lb +. sqr (pop ()) in 
    let input_scaling = sqr (pop ()) in 
    let input_nl_pow = sqr (pop ()) in
    let weight_prms = 
      let ee_h = 0.01 +. sqr (pop ()) in
      let ei_h = 0.01 +. sqr (pop ()) in
      let ie_h = 0.01 +. sqr (pop ()) in
      let ii_h = 0.01 +. sqr (pop ()) in
      let ee_w = pop () in
      let ei_w = pop () in
      let ie_w = pop () in
      let ii_w = pop () in
      let wp (h,w) = { width=w; height=h } in
      { ee = wp (ee_h, ee_w); 
        ei = wp (ei_h, ei_w);
        ie = wp (ie_h, ie_w);
        ii = wp (ii_h, ii_w) } in
    let sigma_eta_prms = 
      let width = pop () in
      let std_e = pop () in
      let std_i = pop () in
      let rho = 0.5 *. (1. +. tanh (pop ())) in
      { width; std_e; std_i; rho } in
    { input_baseline; input_scaling; input_nl_pow; weight_prms; sigma_eta_prms }

  (* pack everything into a big vector of parameters;
     here we need to be careful to do these operations in the same
     order as in [unpack] above *)
  let pack p =
    let v = Vec.create n_prms in
    let c = ref 0 in
    let pop x = incr c; v.{!c} <- x in
    pop (sqrt (p.input_baseline -. input_baseline_lb));
    pop (sqrt p.input_scaling);
    pop (sqrt p.input_nl_pow);
    let pop_height h = 
      if h <= 0.011 then failwith "please raise the initial weight amplitude"
      else pop (sqrt (h -. 0.01)) in
    (* weight heights *)
    pop_height p.weight_prms.ee.height;
    pop_height p.weight_prms.ei.height;
    pop_height p.weight_prms.ie.height;
    pop_height p.weight_prms.ii.height;
    (* weight widths *)
    pop p.weight_prms.ee.width;
    pop p.weight_prms.ei.width;
    pop p.weight_prms.ie.width;
    pop p.weight_prms.ii.width;
    (* sigma_eta *)
    pop p.sigma_eta_prms.width;
    pop p.sigma_eta_prms.std_e;
    pop p.sigma_eta_prms.std_i;
    let rho_p = 
      let z = 2. *. p.sigma_eta_prms.rho -. 1. in
      0.5 *. (log (1. +. z) -. log (1. -. z)) in
    pop rho_p;
    assert (!c = n_prms);
    v

    (*
    1 - input baseline
    2 - input scaling
    3 - input power
    4 - EE height
    5 - EI height
    6 - IE height
    7 - II height
    8 - EE width
    9 - EI width
    10 - IE width
    11 - II width
    12 - noise width
    13 - noise E std
    14 - noise I std
    15 - noise cor
    *)

  (* upper bounds on parameters *)
  let ub = Vec.init n_prms (function 
      | 8 | 9 | 10 | 11 -> PE.(sqrt 2.)
      | 13 | 14 -> 4.
      | _ -> infinity)

  (* -----------------------------------------------------------------------
     @@   Nonlinear moments                                               @@
     ----------------------------------------------------------------------- *)

  type 'a cache = {
    mu: 'a;
    sigma2: 'a;
    sigma: 'a;
    ratio: 'a;
    phi: 'a;
    psi: 'a;
    nu1: 'a; 
  }

  let%diff normal_pdf x = O.(normal_pdf_const *. exp (0.5 *. neg (sqr x)))

  let%diff cache ~mu ~sigma2 =
    let open O in
    let sigma = sqrt sigma2 in
    let ratio = mu / sigma in
    let phi = normal_pdf (module O) ratio in
    let psi = normal_cdf ratio in
    let nu1 = mu * psi + sigma * phi in
    { mu; sigma2; sigma; ratio; phi; psi; nu1 }

  (* I'm hard-coding n=2 (exponent of the power-law) to make use of smart inlining 
     from the compiler *)
  let%diff nu_fun c = O.(k *. (c.mu * c.nu1 + c.sigma2 * c.psi))
  let%diff gamma_fun c = O.(two_k *. c.nu1)

  (* That's the code for threshold-linear, in case we ever need this *)
  (* let%diff nu_fun c = O.(k *. c.nu1)
     let%diff gamma_fun c = O.(k *. c.psi) *)

  (* -----------------------------------------------------------------------
     @@   Moment flows                                                    @@
     ----------------------------------------------------------------------- *)

  let%diff w weight_prms =
    let open O in
    M.init n n (fun i j ->
        let p, s = if i<=m 
          then (if j<=m then weight_prms.ee, 1. else weight_prms.ei, PE.(-1.))
          else (if j<=m then weight_prms.ie, 1. else weight_prms.ii, PE.(-1.)) in
        let i = (pred i) mod m and j = (pred j) mod m in
        let theta_i = PE.(2. *. pi *. float i /. float m) in
        let theta_j = PE.(2. *. pi *. float j /. float m) in
        let cos_diff_minus_one = PE.(cos (theta_i -. theta_j) -. 1.) in
        (s *. p.height) * exp (cos_diff_minus_one *. inv (sqr p.width))) 

  let%diff sigma_eta prms = 
    let open O in
    let var_e = sqr prms.std_e in
    let var_i = sqr prms.std_i in
    let rho = prms.rho in
    M.init n n (fun i j ->
        let xi = if i<=m then var_e else var_i in
        let xj = if j<=m then var_e else var_i in
        let x = sqrt (xi * xj) in
        let x = if ((i<=m) && (j<=m)) || ((i>m) && (j>m)) then x else rho * x in
        let i = (pred i) mod m and j = (pred j) mod m in
        let theta_i = PE.(2. *. pi *. float i /. float m) in
        let theta_j = PE.(2. *. pi *. float j /. float m) in
        let cos_diff_minus_one = PE.(cos (theta_i -. theta_j) -. 1.) in
        x * exp (cos_diff_minus_one *. inv (sqr prms.width))) 
    |> M.add_float_to_diag 0.01

  let%diff j_mat ~w ~gamma =
    let open O in
    let w_eff = M.add_float_to_diag PE.(-1.) (M.scal_cols w gamma) in
    M.scal_rows_float inv_taus w_eff

  (* evolution operators *)

  let%diff dmu_dt ~w ~mu ~nu ~h =
    let open O in
    let z = (h -:|:| mu) +:|:| (w *:||:| nu) in
    V.fmul inv_taus z

  let%diff new_sigma_star ~j_mat ~sigma_eta ~sigma_star ~j_mat_sigma_star =
    let open O in
    let tmp1 = sigma_star  +:||:|| (dt *.:|| j_mat_sigma_star) in
    let tmp2 = M.scal_rows_float inv_taus sigma_eta in
    ((eps1 *.:||  tmp1) +:||:|| (eps2 *.:|| tmp2))


  let%diff new_sigma ~j_mat ~sigma_star ~j_mat_sigma_star ~sigma =
    let open O in
    let b = dt *.:|| (M.scal_cols_float sigma_star inv_taus) in
    let b = b +:||:|| (PE.(dt *. dt) *.:|| j_mat_sigma_star) in
    let prop = M.add_float_to_diag 1. (dt *.:|| j_mat) in
    (prop *:||:|| (sigma *:||:|| M.trans prop)) +:||:|| b +:||:|| M.trans b 

   (* ratio sd/mean regularizer *)
  let%diff reg_sigma_mu_ratio ~mu ~sigma =
    let open O in
    let inv_ratio = V.init m (fun i ->
        let ca = cache (module O) ~mu:(V.get mu PE.(i+m)) ~sigma2:(M.get sigma PE.(i+m) PE.(i+m)) in
        (1.0 /. ca.ratio)) in
    lambda_reg_ratio *. V.sqr_nrm2 inv_ratio

   (* variance regularizer *)
  let%diff reg_var ~sigma =
    let open O in
    let var_diff = V.init m (fun i -> (M.get sigma PE.(i+m) PE.(i+m) - M.get sigma i i)) in
    lambda_reg_var *. V.sqr_nrm2 var_diff


  (* function to simply evolve the moments *)
  let%diff evolve_moments ~mu_init ~sigma_init ~sigma_eta ~w ~h =
    let open O in
    let rec accumulate ~mu ~sigma ~sigma_star t =
      if t=n_time_bins then (mu, sigma, sigma_star)
      else begin
        let sigma_diag = V.init n (fun i -> (M.get sigma i i)) in
        let ca = vinit n (fun i -> cache (module O) ~mu:(V.get mu i) ~sigma2:(V.get sigma_diag i)) in
        let nu = V.init n (fun i -> nu_fun (module O) (vget ca i)) in
        let gamma = V.init n (fun i -> gamma_fun (module O) (vget ca i)) in
        let j_mat = j_mat (module O) ~w ~gamma in
        (* evolve ! *)
        let mu' = mu +:|:| (dt *.:| dmu_dt (module O) ~w ~mu ~nu ~h) in
        let j_mat_sigma_star = j_mat *:||:|| sigma_star in
        let sigma' = new_sigma (module O) ~j_mat ~sigma_star ~j_mat_sigma_star ~sigma in
        let sigma_star' = new_sigma_star (module O) ~j_mat ~sigma_eta ~sigma_star ~j_mat_sigma_star in
        accumulate ~mu:mu' ~sigma:sigma' ~sigma_star:sigma_star' PE.(t+1)
      end in
    (* define sensible initial conditions *)
    let mu0 = match mu_init with Some z -> z | None -> V.of_float (Vec.make n 0.0) in
    let sigma0 = match sigma_init with Some z -> z | None -> M.of_float (F.(4. *.:|| Mat.identity n)) in
    (* work out the sigma_star that's consistent with the current sigma *)
    (* mmmh... for now, pick a dump sigma_star *)
    let sigma_star0 = M.of_float (F.(4. *.:|| Mat.identity n)) in
    (mu0, sigma0, sigma_star0), accumulate ~mu:mu0 ~sigma:sigma0 ~sigma_star:sigma_star0 0

  let step_fun i = if i<(n_time_bins - min_time_bins) then 0. else 1. 
  let temporal_weighting = Vec.init n_time_bins (fun i -> step_fun i)

  let%diff evolution_costs ~sigma_eta ~w ~h ~target_mu ~target_sigma = 
    let open O in
    let rec accumulate ~mu ~sigma ~sigma_star ((accu_mean, accu_var, accu_cov) as accu) t =
      if t =n_time_bins then (accu, mu, sigma, sigma_star, reg_sigma_mu_ratio (module O) ~mu ~sigma,  reg_var (module O) ~sigma)
      else begin
        let sigma_diag = V.init n (fun i -> (M.get sigma i i)) in
        let ca = vinit n (fun i -> cache (module O) ~mu:(V.get mu i) ~sigma2:(V.get sigma_diag i)) in
        let nu = V.init n (fun i -> nu_fun (module O) (vget ca i)) in
        let gamma = V.init n (fun i -> gamma_fun (module O) (vget ca i)) in
        let j_mat = j_mat (module O) ~w ~gamma in
        (* evolve ! *)
        let mu' = mu +:|:| (dt *.:| dmu_dt (module O) ~w ~mu ~nu ~h) in
        let j_mat_sigma_star = j_mat *:||:|| sigma_star in
        let sigma' = new_sigma (module O) ~j_mat ~sigma_star ~j_mat_sigma_star ~sigma in
        let sigma_star' = new_sigma_star (module O) ~j_mat ~sigma_eta ~sigma_star ~j_mat_sigma_star in
        (* compute the distance from target *)
        let accu' = if t mod subsamp_bins = 0 then (
            let wt = V.getf temporal_weighting (succ t) in
            let cost_mean = lambda_mean *. V.sqr_nrm2 (V.init m (fun i -> target_mu.{i} -. V.get mu i)) in
            let cost_var = lambda_var *. V.sqr_nrm2 (V.init m (fun i -> target_sigma.{i,i} -. M.get sigma i i)) in
            let cost_cov = lambda_cov *. M.sqr_frob (M.init m m (fun i j ->
                target_sigma.{i,j} -. M.get sigma i j)) in
            accu_mean + wt * cost_mean, accu_var + wt * cost_var, accu_cov + wt * cost_cov
          ) else accu in
        accumulate ~mu:mu' ~sigma:sigma' ~sigma_star:sigma_star' accu' PE.(t+1)
      end in
    (* define sensible initial conditions *)
    let accu = of_float 0., of_float 0., of_float 0. in
    let mu = V.of_float (Vec.make n 0.0) in
    let sigma = M.of_float (F.(4. *.:|| Mat.identity n)) in
    let sigma_star = M.of_float (F.(4. *.:|| Mat.identity n)) in
    accumulate ~mu ~sigma ~sigma_star accu 0

  let%diff h prms h_vec = 
    let open O in
    V.init n (fun i -> 
        prms.input_scaling * exp (prms.input_nl_pow * log (h_vec.{i} +. prms.input_baseline)))

  let%diff reg_sigma_eta sigma_eta =
    let open O in
    lambda_noise *. V.sqr_nrm2 (V.init n (fun i -> input_noise_size -. M.get sigma_eta i i))

  let%diff reg_sigma_eta_width prms =
    let open O in
    lambda_noise_width *. inv (sqr prms.width)

  let%diff reg_rates mu =
    let open O in
    let z = ref (of_float 0.) in
    for i=1 to m  do z := !z + (V.get mu i) done;
    for i=succ m to n do z := !z - (V.get mu i) done;
    lambda_rates *. sqr !z

  let%diff reg_sigma sigma =
    let open O in
    let z = ref (of_float 0.) in
    for i=1 to m do z := !z + (M.get sigma i i) done;
    for i=succ m to n do z := !z - (M.get sigma i i) done;
    lambda_sigma *. sqr !z

  (* slowness cost *)
  let%diff slowness ~w ~mu ~sigma =
    let open O in
    let gamma = V.init n (fun i ->
        let ca = cache (module O) ~mu:(V.get mu i) ~sigma2:(M.get sigma i i) in
        gamma_fun (module O) ca) in
    let j_mat = j_mat (module O) ~w ~gamma in
    let rec iterate t cost s = 
      if t=n_time_bins_slow then lambda_slow *. cost
      else begin
        let s_norm_diag = V.init n (fun i -> M.get s i i / M.get sigma i i) in 
        let new_cost = cost + V.sqr_nrm2 s_norm_diag in (* + V.sum s_norm_diag in *)
        let new_s = s +:||:|| (dt *.:|| (s *:||:|| M.trans j_mat)) in
        iterate (succ t) new_cost new_s
      end in
    iterate 0 (of_float 0.) sigma

  let%diff objectives p =
    let open O in
    let prms = unpack (module O) p in
    let w = w (module O) prms.weight_prms in
    let sigma_eta = sigma_eta (module O) prms.sigma_eta_prms in
    let noise_reg = reg_sigma_eta (module O) sigma_eta in
    let noise_width_reg = reg_sigma_eta_width (module O) prms.sigma_eta_prms in
    (* looping through all targets *)
    noise_reg, noise_width_reg, Array.map (fun targ ->
        let h = h (module O) prms targ.h_vec in
        let ((cost_mean, cost_var, cost_cov), mu, sigma, sigma_star,reg_ratio,reg_variance) = evolution_costs (module O) ~sigma_eta ~w ~h ~target_mu:targ.mu_vec ~target_sigma:targ.sigma_mat in 
        let cost_slow = slowness (module O) ~w ~mu ~sigma in
        ((cost_mean, cost_var, cost_cov, cost_slow, reg_ratio, reg_variance), mu, sigma, sigma_star)
      ) targets

  let%diff objective p =
    let open O in
    let noise_reg, noise_width_reg, costs = objectives (module O) p in
    let evolution_cost = Array.fold_left (fun accu ((cost_mean, cost_var, cost_cov, cost_slow, reg_ratio,reg_variance), mu, sigma, _) -> 
        accu + cost_mean + cost_var + cost_cov + cost_slow + reg_ratio + reg_variance) (of_float 0.) costs in
    evolution_cost + noise_reg + noise_width_reg

  let%diff objective_slowness_only p =
    let open O in
    let noise_reg, noise_width_reg, costs = objectives (module O) p in
    Array.fold_left (fun accu ((_, _, _, cost_slow,reg_ratio,reg_variance), _, _, _) -> 
        accu + cost_slow) (of_float 0.) costs


  (* ---------- sample-based approach ----------- *)

  module Samples = struct

    type t_noise_prms = { 
      n_trials: int;
      noise_bank: (float, prec, Bigarray.fortran_layout) Bigarray.Array3.t;
      eta_init: mat;
      ones: vec;
    }

    let noise_prms = match n_trials with
      | None -> None
      | Some n_trials -> 
        (* noise matrix *)
        let noise_bank = Bigarray.(Array3.create Float64 fortran_layout n n_trials n_time_bins) in
        let eta_init =  Mat.make0 n n_trials in
        let ones = Vec.make n_trials 1. in
        Some { n_trials; noise_bank; eta_init; ones }

    let redraw_noise () = match noise_prms with
      | None -> ()
      | Some np ->
        let open Bigarray in
        for t=1 to n_time_bins do
          for k=1 to np.n_trials do
            let v = Array3.slice_right_1 np.noise_bank k t in
            for i=1 to n do 
              v.{i} <- Stuff.Rand.gaussian_noise 1.0;
              np.eta_init.{i,k} <- Stuff.Rand.gaussian_noise 1.0
            done
          done
        done


    (* ell = cholesky factor of sigma_eta *)
    let%diff eta ~noise_prms ell =
      let open O in
      let x = ref (ell *:||.|| noise_prms.eta_init) in
      Array.init n_time_bins (fun t ->  (* maybe this could be sped up by stacking the time frames
                                           and using a single matrix-matrix product (?) *)
          x := (eps1 *.:|| !x) +:||:|| (eps3 *.:|| (ell *:||.|| Bigarray.Array3.slice_right_2 noise_prms.noise_bank (succ t)));
          !x)


    let temporal_weighting = Vec.make n_time_bins 1.

    let%diff evolution_costs ~noise_prms ~eta ~w ~h ~target_mu ~target_sigma =
      let open O in
      let rec accumulate u ((accu_mean, accu_var, accu_cov) as accu) t =
        (* update u *)
        let r = M.map_threshold_quadratic ~k:P.k u in
        let du = (M.add_col_vec (M.neg u) h) +:||:|| (w *:||:|| r) +:||:|| eta.(t) in
        let u = u +:||:|| (M.scal_rows_float dt_inv_taus du) in
        (* soft threshold u *)
        let soft_gain = 100. in
        let u = soft_gain *.:|| M.map_tanh (PE.(1. /. soft_gain) *.:|| u) in
        (* compute the sample moment estimates *)
        let mu = PE.(1.0 /. float noise_prms.n_trials) *.:| (u *:||.| noise_prms.ones) in
        let u_centered = M.add_col_vec u (V.neg mu) in
        let sigma = PE.(1.0 /. float noise_prms.n_trials) *.:|| (u_centered *:||:|| M.trans u_centered) in
        let accu = if t mod subsamp_bins = 0 then (
            let wt = V.getf temporal_weighting (succ t) in
            let cost_mean = lambda_mean *. V.sqr_nrm2 (V.init m (fun i -> target_mu.{i} -. V.get mu i)) in
            let cost_var = lambda_var *. V.sqr_nrm2 (V.init m (fun i -> target_sigma.{i,i} -. M.get sigma i i)) in
            let cost_cov = lambda_cov *. M.sqr_frob (M.init m m (fun i j -> target_sigma.{i,j} -. M.get sigma i j)) in
            accu_mean + wt * cost_mean, accu_var + wt * cost_var, accu_cov + wt * cost_cov
          ) else accu in
        if t = pred n_time_bins then (accu, (mu, sigma))
        else accumulate u accu (succ t) 
      in
      let u0 = M.of_float (Mat.init n noise_prms.n_trials (fun _ _ -> Stuff.Rand.gaussian_noise 2.)) in
      let accu = of_float 0., of_float 0., of_float 0. in
      accumulate u0 accu 0


    let%diff objectives ~noise_prms p =
      let open O in
      let prms = unpack (module O) p in
      let w = w (module O) prms.weight_prms in
      let sigma_eta = sigma_eta (module O) prms.sigma_eta_prms in
      let noise_width_reg = reg_sigma_eta_width (module O) prms.sigma_eta_prms in
      let noise_reg = reg_sigma_eta (module O) sigma_eta in
      let ell = M.chol sigma_eta in
      let eta = eta (module O) ~noise_prms ell in
      (* looping through all targets *)
      noise_reg, noise_width_reg, Array.map (fun targ ->
          let h = h (module O) prms targ.h_vec in
          evolution_costs (module O) ~noise_prms ~eta 
            ~w ~h ~target_mu:targ.mu_vec ~target_sigma:targ.sigma_mat
        ) targets

    let%diff objective p =
      let open O in
      match noise_prms with
      | None -> assert false
      | Some noise_prms ->
        let noise_reg, noise_width_reg, costs = objectives (module O) ~noise_prms p in
        let evolution_cost = Array.fold_left (fun accu ((cost_mean, cost_var, cost_cov), _) -> 
            accu + cost_mean + cost_var + cost_cov) 
            (of_float 0.) costs in
        evolution_cost + noise_reg + noise_width_reg

  end

  let f = { n=n_prms; f=(match P.n_trials with Some _ -> Samples.objective | None -> objective)}
  let f_slowness = { n=n_prms; f=objective_slowness_only }

end

(* -----------------------------------------------------------------------
   @@   Some convenience functions                                      @@
   ----------------------------------------------------------------------- *)


(* generate a random covariance matrix and the associated upper-tri Cholesky factor *)
let random_sigma_ell ~avg_var ~cor_dist_std n =
  let sigma = Stuff.Wishart.cov_mat_inv_wishart ~avg_var ~cor_dist_std n
              |> Mat.of_array in
  let ell = Mat.chol `Upper sigma in
  sigma, ell

