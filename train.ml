open Printf
open Math.Core
open Autodiff
open Objective
module F = Operators.On_floats

let () = Util.set_verbose false


let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [directory]")
let m = Cmdargs.(get_int "-m" |> force ~usage:"-m [number of ring sites]")
let n_targ = Cmdargs.(get_int "-n_targ" |> force ~usage:"-n [number of targets]")
let save_every = Cmdargs.(get_int "-save_every" |> default 10)
let max_iter = Cmdargs.(get_int "-max_iter" |> default 1000)
let rand_init = Cmdargs.(get_bool "-rand_init" |> default false)


let in_dir = sprintf "%s/%s" dir

(* read the targets *)
let read_targets () = Array.init n_targ (fun k ->
    printf "k = %i\n%!" k;
    let h_vec = Vec.read (in_dir (sprintf "h%i" k)) (2*m) in
    let mu_vec = Vec.read (in_dir (sprintf "mu%i" k)) m in
    let sigma_mat = Mat.read (in_dir (sprintf "sigma%i" k)) m m in
    { h_vec; mu_vec; sigma_mat })

module P = struct
  let verbose = true
  let m = m
  let k = 0.3
  let tau_e = Vec.make m 20E-3
  let tau_i = Vec.make m 10E-3
  let tau_eta = Cmdargs.(get_float "-tau_eta" |> default 50E-3)
  let input_noise_size = Cmdargs.(get_float "-input_noise_size" |> default 4.0)
  let lambda_mean = Cmdargs.(get_float "-lambda_mean" |> default 1.)
  let lambda_var = Cmdargs.(get_float "-lambda_var" |> default 1.)
  let lambda_cov = Cmdargs.(get_float "-lambda_cov" |> default 1.)
  let lambda_slow = Cmdargs.(get_float "-lambda_slow" |> default 0.)
  let lambda_noise = Cmdargs.(get_float "-lambda_noise" |> default 0.)
  let lambda_noise_width = Cmdargs.(get_float "-lambda_noise_width" |> default 0.)
  let lambda_rates = Cmdargs.(get_float "-lambda_rates" |> default 0.0)
  let lambda_sigma = Cmdargs.(get_float "-lambda_sigma" |> default 0.0)
  let lambda_reg_ratio = Cmdargs.(get_float "-lambda_reg_ratio" |> default 0.0)
  let lambda_reg_var = Cmdargs.(get_float "-lambda_reg_var" |> default 0.0)
  let targets = read_targets ()
  let input_baseline_lb =
    let m = targets |> Array.fold_left (fun accu targ -> min accu (Vec.min targ.h_vec)) 0. in
    0.1 -. m
  let dt = Cmdargs.(get_float "-dt" |> default 0.2E-3)
  let t_max = Cmdargs.(get_float "-t_max" |> default 0.1)
  let t_max_slow = Cmdargs.(get_float "-t_max_slow" |> default 0.1)
  let t_subsamp = Cmdargs.(get_float "-t_subsamp" |> default 10.0E-3)
  let n_time_bins_count = Cmdargs.(get_int "-n_time_bins_count" |> default 10)
  let n_trials = Cmdargs.get_int "-n_trials"
end

(* used for optimization *)
module X = Make (P)

(* used for periodically testing the forward evolution of the moments *)
module XTest = Make (struct include P
    let verbose = false
    let t_max = 1.
    let dt = 0.2E-3
  end)
open X

let () = X.Samples.redraw_noise ()


(* pick a sensible initial parameter vector *)
let x = match Cmdargs.get_string "-reuse" with
  | Some file_name ->
    printf "I am REUSING a previously saved state vector\n%!";
    Vec.read_bin file_name
  | None ->
    if rand_init then
      let () = printf "I am GENERATING a random initial state vector.\n%!" in
      let () = Random.self_init() in
      let input_baseline = P.input_baseline_lb +. Random.float 2. in
      let input_scaling = (0.5 +. Random.float 1.) in
      let input_nl_pow = 1. in
      let weight_prms = {
        ee = X.{ width=0.8; height=0.02 +. Random.float 0.02 };
        ei = X.{ width=0.8; height=0.02 +. Random.float 0.02 };
        ie = X.{ width=0.8; height=0.02 +. Random.float 0.02 };
        ii = X.{ width=0.8; height=0.02 +. Random.float 0.02 } } in
      let sigma_eta_prms = { width=0.8; std_e=2.; std_i=2.; rho = 0.8 } in
      pack { input_baseline; input_scaling; input_nl_pow; weight_prms; sigma_eta_prms }
    else
      let () = printf "I am GENERATING a sensible predetermined initial state vector.\n%!" in
      let input_baseline = P.input_baseline_lb +. 1. in
      let input_scaling = 1. in
      let input_nl_pow = 1. in
      let weight_prms = {
        ee = X.{ width=0.8; height=0.02 };
        ei = X.{ width=0.8; height=0.02 };
        ie = X.{ width=0.8; height=0.02 };
        ii = X.{ width=0.8; height=0.02 } } in
      let sigma_eta_prms = { width=0.8; std_e=2.; std_i=2.; rho = 0.8 } in
      (* pack it up *)
      pack { input_baseline; input_scaling; input_nl_pow; weight_prms; sigma_eta_prms }


let save_results x label =
  Vec.save_bin (in_dir "state_vector.bin") x;
  let in_dir s = in_dir (sprintf "%s_%s" s label) in
  let prms = X.unpack (module F) x in
  let w = w (module F) prms.weight_prms in
  let w_ee = prms.weight_prms.ee in
  let w_ei = prms.weight_prms.ei in
  let w_ii = prms.weight_prms.ii in
  let w_ie = prms.weight_prms.ie in
  let sigma_eta = sigma_eta (module F) prms.sigma_eta_prms in
  let var_width = prms.sigma_eta_prms.width in
  let var_e = sqr prms.sigma_eta_prms.std_e in
  let var_i = sqr prms.sigma_eta_prms.std_i in
  let rho = prms.sigma_eta_prms.rho in
  Mat.save (in_dir "w") w;
  Vec.save (in_dir "w_ee_width") (Vec.make 1 w_ee.width);
  Vec.save (in_dir "w_ei_width") (Vec.make 1 w_ei.width);
  Vec.save (in_dir "w_ii_width") (Vec.make 1 w_ii.width);
  Vec.save (in_dir "w_ie_width") (Vec.make 1 w_ie.width);
  Vec.save (in_dir "w_ee_height") (Vec.make 1 w_ee.height);
  Vec.save (in_dir "w_ei_height") (Vec.make 1 w_ei.height);
  Vec.save (in_dir "w_ii_height") (Vec.make 1 w_ii.height);
  Vec.save (in_dir "w_ie_height") (Vec.make 1 w_ie.height);
  Mat.save (in_dir "sigma_eta") sigma_eta;
  Vec.save (in_dir "var_width") (Vec.make 1 var_width);
  Vec.save (in_dir "var_e") (Vec.make 1 var_e);
  Vec.save (in_dir "var_i") (Vec.make 1 var_i);
  Vec.save (in_dir "rho") (Vec.make 1 rho);
  Vec.save (in_dir "input_scaling") (Vec.make 1 prms.input_scaling);
  Vec.save (in_dir "input_nl_pow") (Vec.make 1 prms.input_nl_pow);
  Vec.save (in_dir "input_baseline") (Vec.make 1 prms.input_baseline);
  Array.iteri (fun k tk ->
      let h = h (module F) prms tk.h_vec in
      Vec.save (in_dir (sprintf "h_true_%i" k)) h;
    ) P.targets


let test_forward_evolution x label =
  let open XTest in
  let in_dir s = in_dir (sprintf "%s_%s" s label) in
  let prms = unpack (module F) x in
  let w = w (module F) prms.weight_prms in
  let sigma_eta = sigma_eta (module F) prms.sigma_eta_prms in
  (* looping through all targets *)
  Array.iteri (fun k targ ->
      let h = h (module F) prms targ.h_vec in
      let costs, mu, sigma, _, e, v = evolution_costs (module F)
          ~w ~h ~sigma_eta ~target_mu:targ.mu_vec ~target_sigma:targ.sigma_mat in
      Vec.save (in_dir (sprintf "mu_%i" k)) mu;
      Mat.save (in_dir (sprintf "sigma_%i" k)) sigma;
      Vec.save (in_dir (sprintf "std_%i" k)) (Vec.map sqrt (Mat.diag sigma));
    ) P.targets


let () = 
  save_results x "init"
(* test_forward_evolution x "evolved_init" *)

(* optimize ! *)

(* gradient test *)
let () =  if Cmdargs.check "-test_gradients" then begin
    printf "Testing gradients...%!";
    let open X in
    assert (test_gradient f `Reverse x);
    printf " done.\n%!"
  end


(* use the lbfgs library to minimize our cost function *)
let _, fdf = grad X.f
let f_df x g =
  let cost, g_ = fdf x in
  Vec.blit g_ g;
  cost

(* gradient of the slowness cost, from which we remove the component aligned onto
   the gradient of the objective only (without slowness) *)
let f_df_slowness_perp =
  let _, fdf_slowness = grad X.f_slowness in
  fun x g ->
    let _, g_total = fdf x in
    let cost, g_slowness = fdf_slowness x in
    let g_obj = F.(g_total -:|:| g_slowness) in
    let g_obj = Vec.normalize 1. g_obj in
    Vec.blit F.(g_slowness -:|:| ((dot g_slowness g_obj) *.:| g_obj)) g;
    cost


let common_stop k cost =
  if k mod 10 = 0 then Gc.compact (); (* to avoid running into swap memory *)
  Vec.save ~append:true (in_dir "running_cost") (Vec.of_array [| cost |]);
  (* save more things every 10 iterations *)
  if k mod save_every = 0 then begin

    (* save the ADF moments *)
    let noise_reg, noise_width_reg, data = X.objectives (module F) x in
    (* total cost *)
    let total_cost = data |> Array.map (fun ((x,y,z,s,r,v), _, _, _) -> x +. y +. z +. s +. r +. v) |> Stuff.Core.sum in
    let total_slowness_cost = data |> Array.map (fun ((_,_,_,s,_,_), _, _, _) -> s) |> Stuff.Core.sum in
    (* pattern-specific sub-costs *)
    let all_costs = data |> Array.map (fun ((x,y,z,s,r,v), _, _, _) -> [| x; y; z; s; r; v |])
                    |> Stuff.Core.trans |> Stuff.Core.flatten in
    Mat.save ~append:true (in_dir "costs")
      (Mat.of_array ([| Array.append [| total_cost; total_slowness_cost; noise_reg; noise_width_reg |] all_costs |]));
    (* save the final moments *)
    Array.iteri (fun k (_, mu, sigma, sigma_star) ->
        Vec.save (in_dir (sprintf "mu_learn_%i" k)) mu;
        Mat.save (in_dir (sprintf "sigma_learn_%i" k)) sigma;
        Vec.save (in_dir (sprintf "sigma_diag_learn_%i" k)) (Mat.diag sigma);
        Mat.save (in_dir (sprintf "sigma_star_learn_%i" k)) sigma_star
      ) data;

    (* save the sample-based moments *)
    begin match X.Samples.noise_prms with
      | None -> ()
      | Some noise_prms ->
        let noise_reg, noise_width_reg, data = X.Samples.objectives (module F) ~noise_prms x in
        Array.iteri (fun k (_, (mu, sigma)) ->
            Vec.save (in_dir (sprintf "mu_learn_samples_%i" k)) mu;
            Mat.save (in_dir (sprintf "sigma_learn_samples_%i" k)) sigma;
            Vec.save (in_dir (sprintf "sigma_diag_learn_samples_%i" k)) (Mat.diag sigma);
          ) data;
        Vec.save (in_dir "temporal_weighting") X.Samples.temporal_weighting;
    end;
  end;
  if k mod save_every = 0 then save_results x "learn";
  if k mod save_every = 0 then test_forward_evolution x "evolved";
  (* we'll stop manually *)
  if (k >= max_iter) then true else false


(* change the weighting function used in weighting different time points during the transient;
   if -no_progression is provided, will go straight to the final weighting *)
let update_temporal_weighting iter =
  let zero_up_to = min (X.n_time_bins - X.min_time_bins) 
      (if Cmdargs.check "-no_progression" then max_int else round (float n_time_bins *. float iter /. 200.)) in
  for i=1 to X.n_time_bins do X.Samples.temporal_weighting.{i} <- if i<zero_up_to then 0. else 1. done



let () = match P.n_trials with

  | None -> (* BFGS *)
    printf "MOMENT BASED OPTIMIZATION\n\n\n%!";
    if Cmdargs.check "-slowness_only" then begin

      let stop k cost = 
        printf "iteration %5i | cost = %.5f\n%!" k cost;
        common_stop k cost in
      let eta = Cmdargs.(get_float "-eta" |> default 0.002) in
      let epsilon = Cmdargs.(get_float "-epsilon" |> default 10E-8) in
      let beta1 = Cmdargs.(get_float "-beta1" |> default 0.9) in
      let beta2 = Cmdargs.(get_float "-beta2" |> default 0.999) in
      let clip = sqrt (float X.n_prms) *. 5. in
      Adam.min ~eta ~epsilon ~beta1 ~beta2 ~ub:X.ub ~clip ~stop f_df_slowness_perp x |> ignore

    end else begin
      let stop st =
        let k = Lbfgs.iter st in
        let cost = Lbfgs.previous_f st in
        common_stop k cost in
      let rec attempt () =
        try Lbfgs.(F.min ~print:(Every 1) ~u:X.ub ~factr:1E1 ~pgtol:0. ~corrections:20 ~stop f_df x |> ignore)
        with _ ->  begin
            save_results x "failed";
            (* take a mini gradient step *)
            let _, g = fdf x in
            Vec.sub ~z:x x F.(0.001 *.:| g) |> ignore;
            attempt ()
          end
      in attempt ()
    end

  | Some _ -> (* ADAM *)
    update_temporal_weighting 1;
    printf "SAMPLE BASED OPTIMIZATION\n\n\n%!";
    let stop k cost = 
      printf "iteration %5i | cost = %.5f\n%!" k cost;
      X.Samples.redraw_noise ();
      update_temporal_weighting k;
      common_stop k cost in
    let eta = Cmdargs.(get_float "-eta" |> default 0.002) in
    let epsilon = Cmdargs.(get_float "-epsilon" |> default 10E-8) in
    let beta1 = Cmdargs.(get_float "-beta1" |> default 0.9) in
    let beta2 = Cmdargs.(get_float "-beta2" |> default 0.999) in
    let clip = sqrt (float X.n_prms) *. 5. in
    Adam.min ~eta ~epsilon ~beta1 ~beta2 ~ub:X.ub ~clip ~stop f_df x |> ignore

let () = save_results x "final"; test_forward_evolution x "evolved_final"

