"""TubulePy.py

GUI to compute protofilament skew angle theta for any N_S microtubule configuration
using the Lattice Accommodation Model (LAM) formula:

    theta = arctan(1/dx * (S * a / N) - r)

If used in a publication, please cite:
Chrétien and Fuller, J. Mol. Biol., 2000, 298:663-676.

The GUI accepts N, S, a, r, dx (single or ranges) and shows results.
"""
import tkinter as tk
from tkinter import messagebox, ttk, filedialog
import math
import os
import subprocess
import time
import json
import sys
from pathlib import Path

# Store Tk image references to prevent garbage collection
png_thumbnails = []
from datetime import datetime

APP_NAME = "TubulePy"
CHIMERAX_ENV_VAR = "TUBULEPY_CHIMERAX"
CHIMERAX_DEFAULT_BIN = "/Applications/ChimeraX-1.11.app/Contents/bin/ChimeraX"


def _config_file():
    base = None
    if os.name == "nt":
        base = os.getenv("APPDATA")
    elif sys.platform == "darwin":
        base = os.path.join(Path.home(), "Library", "Application Support")
    else:
        base = os.path.join(Path.home(), ".config")
    return Path(base) / APP_NAME / "config.json"


def load_chimerax_bin():
    env_path = os.getenv(CHIMERAX_ENV_VAR)
    if env_path and os.path.exists(env_path):
        return env_path
    cfg = _config_file()
    if cfg.exists():
        try:
            data = json.loads(cfg.read_text())
            saved = data.get("chimerax_bin")
            if saved and os.path.exists(saved):
                return saved
        except Exception:
            pass
    return CHIMERAX_DEFAULT_BIN


def save_chimerax_bin(path_str):
    try:
        cfg = _config_file()
        cfg.parent.mkdir(parents=True, exist_ok=True)
        cfg.write_text(json.dumps({"chimerax_bin": path_str, "saved_at": datetime.now().isoformat()}))
    except Exception:
        pass


CHIMERAX_BIN = load_chimerax_bin()

# Widgets that need broader scope
generated_combo = None
latest_chimerax_run_dir = ""

# Try to import matplotlib, but make it optional
try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # noqa: F401 (not always used)
    MATPLOTLIB_AVAILABLE = True
except ImportError as e:
    MATPLOTLIB_AVAILABLE = False
    MATPLOTLIB_ERROR = str(e)

# Optional dependencies for MRC projection
try:
    import mrcfile
    MRCFILE_AVAILABLE = True
except ImportError as e:
    mrcfile = None
    MRCFILE_AVAILABLE = False
    MRCFILE_ERROR = str(e)

try:
    import numpy as np
    NP_AVAILABLE = True
except ImportError as e:
    np = None
    NP_AVAILABLE = False
    NP_ERROR = str(e)

try:
    from PIL import Image, ImageTk
    PIL_AVAILABLE = True
except ImportError:
    Image = None
    ImageTk = None
    PIL_AVAILABLE = False


def compute_theta_single():
    """Compute theta for single N, S values."""
    try:
        N = int(var_N.get())
        S = float(var_S.get())
        a = float(var_a.get())
        r = var_r_exact  # use exact stored value
        dx = float(var_dx.get())
        a_for_outputs = a * 2 if var_use_dimer.get() else a

        if N == 0:
            raise ValueError("N must be non-zero")
        if dx == 0:
            raise ValueError("dx must be non-zero")

        inner = (1.0 / dx) * ((S * a / N) - r)
        theta_rad = math.atan(inner)
        theta_deg = math.degrees(theta_rad)

        lbl_result.config(text=f"θ = {theta_deg:.6f}°")

        n_s_config = f"{N}_{S:.1f}"
        theta_rad_for_formulas = math.radians(theta_deg)
        r_su_A = math.cos(theta_rad_for_formulas) * a_for_outputs
        Phi_su_deg = (math.sin(2 * theta_rad_for_formulas) * a_for_outputs * 180) / (N * dx)
        Delta_C_A = (dx / (2 * math.pi)) * ((N / math.cos(theta_rad_for_formulas)) - 13)
        r_pf_A = math.sin(theta_rad_for_formulas) * dx + r * math.cos(theta_rad_for_formulas)
        if var_use_dimer.get():
            Phi_pf_deg = (S * Phi_su_deg - 720) / (2 * N)
        else:
            Phi_pf_deg = (S * Phi_su_deg - 360) / N

        try:
            if not var_add_refs.get():
                for it in table.get_children():
                    table.delete(it)
            table.insert(
                "",
                "end",
                text="☑",
                values=(
                    n_s_config,
                    N,
                    S,
                    f"{theta_deg:.4f}",
                    f"{r_su_A:.4f}",
                    f"{Phi_su_deg:.4f}",
                    f"{Delta_C_A:.4f}",
                    f"{r_pf_A:.4f}",
                    f"{Phi_pf_deg:.4f}",
                ),
            )
        except Exception:
            pass

    except ValueError as e:
        messagebox.showerror("Input error", f"Invalid input: {e}")
    except Exception as e:
        messagebox.showerror("Error", str(e))


def on_r_focus_out(event):
    """When the r entry loses focus, parse input and store exact value."""
    global var_r_exact
    try:
        val = float(var_r.get())
        var_r_exact = val
        var_r.set(f"{var_r_exact:.2f}")
    except Exception:
        messagebox.showerror("Input error", "Invalid value for r. Reverting to previous value.")
        var_r.set(f"{var_r_exact:.2f}")


def compute_theta_range():
    """Compute theta over ranges of N and S, filtering by theta range (in degrees)."""
    try:
        N_min = int(var_N_min.get())
        N_max = int(var_N_max.get())
        S_min = float(var_S_min.get())
        S_max = float(var_S_max.get())
        a = float(var_a.get())
        r = var_r_exact
        dx = float(var_dx.get())
        a_for_outputs = a * 2 if var_use_dimer.get() else a
        theta_min_deg = float(var_theta_min.get())
        theta_max_deg = float(var_theta_max.get())

        theta_min_rad = math.radians(theta_min_deg)
        theta_max_rad = math.radians(theta_max_deg)

        if N_min <= 0 or N_max <= 0:
            raise ValueError("N values must be > 0")
        if S_min <= 0 or S_max <= 0:
            raise ValueError("S values must be > 0")
        if dx == 0:
            raise ValueError("dx must be non-zero")

        for item in table.get_children():
            table.delete(item)

        count = 0
        for N in range(N_min, N_max + 1):
            S = S_min
            while S <= S_max:
                inner = (1.0 / dx) * ((S * a / N) - r)
                theta_rad = math.atan(inner)
                theta_deg = math.degrees(theta_rad)

                should_insert = True
                if var_use_theta_filter.get():
                    should_insert = theta_min_rad <= theta_rad <= theta_max_rad

                if should_insert:
                    n_s_config = f"{N}_{S:.1f}"
                    theta_rad_for_formulas = math.radians(theta_deg)
                    r_su_A = math.cos(theta_rad_for_formulas) * a_for_outputs
                    Phi_su_deg = (math.sin(2 * theta_rad_for_formulas) * a_for_outputs * 180) / (N * dx)
                    Delta_C_A = (dx / (2 * math.pi)) * ((N / math.cos(theta_rad_for_formulas)) - 13)
                    r_pf_A = math.sin(theta_rad_for_formulas) * dx + r * math.cos(theta_rad_for_formulas)
                    if var_use_dimer.get():
                        Phi_pf_deg = (S * Phi_su_deg - 720) / (2 * N)
                    else:
                        Phi_pf_deg = (S * Phi_su_deg - 360) / N

                    table.insert(
                        "",
                        "end",
                        text="☑",
                        values=(
                            n_s_config,
                            N,
                            S,
                            f"{theta_deg:.4f}",
                            f"{r_su_A:.4f}",
                            f"{Phi_su_deg:.4f}",
                            f"{Delta_C_A:.4f}",
                            f"{r_pf_A:.4f}",
                            f"{Phi_pf_deg:.4f}",
                        ),
                    )
                    count += 1
                S += 0.5

        filter_status = "(filtered by θ range)" if var_use_theta_filter.get() else "(no filter)"
        lbl_result.config(text=f"Computed {count} values {filter_status}")

    except ValueError as e:
        messagebox.showerror("Input error", f"Invalid input: {e}")
    except Exception as e:
        messagebox.showerror("Error", str(e))


def plot_n_vs_theta():
    """Plot N vs theta from the current table data."""
    if not MATPLOTLIB_AVAILABLE:
        error_msg = (
            "Matplotlib is not available in the current Python environment.\n\n"
            f"Error details: {MATPLOTLIB_ERROR}\n\n"
            "To fix this, ensure you're using the correct Python interpreter:\n"
            "1. Open a terminal\n"
            "2. Run: source /Users/Denis_1/Documents/VS_Code/.venv/bin/activate\n"
            "3. Run: python /Users/Denis_1/Documents/VS_Code/LAM_Simple/LamSimple.py\n\n"
            "Or set VS Code's Python interpreter to the venv:\n"
            "Cmd+Shift+P → Python: Select Interpreter → Choose .venv/bin/python"
        )
        messagebox.showerror("Matplotlib not available", error_msg)
        return

    try:
        items = table.get_children()
        n_values_blue, theta_values_blue = [], []
        n_values_red, theta_values_red = [], []

        for item in items:
            keep_text = table.item(item, "text")
            if str(keep_text) not in ("☑", "✓", "True", "1"):
                continue
            values = table.item(item)["values"]
            try:
                n_val = float(values[1])
            except Exception:
                try:
                    n_val = float(str(values[0]).split("_")[0])
                except Exception:
                    continue
            try:
                theta_val = float(values[3])
            except Exception:
                try:
                    theta_val = float(str(values[3]))
                except Exception:
                    continue
            try:
                s_val = float(values[2])
                s_twice = round(s_val * 2)
                is_half = (abs(s_val * 2 - s_twice) < 1e-6) and (s_twice % 2 == 1)
            except Exception:
                is_half = False

            if is_half:
                n_values_red.append(n_val)
                theta_values_red.append(theta_val)
            else:
                n_values_blue.append(n_val)
                theta_values_blue.append(theta_val)

        fig, ax = plt.subplots(figsize=(8, 6))
        if n_values_blue:
            ax.scatter(n_values_blue, theta_values_blue, s=50, alpha=0.6, color="blue", label="S ends with .0")
        if n_values_red:
            ax.scatter(n_values_red, theta_values_red, s=50, alpha=0.6, color="red", label="S ends with .5")
        ax.set_xlabel("N (protofilaments)", fontsize=12)
        ax.set_ylabel("θ (degrees)", fontsize=12)
        ax.set_title("Protofilament Number vs Skew Angle", fontsize=14)
        ax.grid(True, alpha=0.3)
        if n_values_blue or n_values_red:
            ax.legend()

        plt.tight_layout()
        plt.show()
    except Exception as e:
        messagebox.showerror("Plot error", f"Error creating plot: {e}")


def export_to_txt():
    """Export table results to a tab-separated TXT file."""
    try:
        items = table.get_children()
        if not items:
            messagebox.showwarning("No data", "Please compute results first.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"Ref_{timestamp}.txt"
        filepath = os.path.expanduser(f"{filename}")

        columns = table["columns"]
        col_alias = {"N_S": "MTtype", "θ (°)": "Theta_deg", "theta (°)": "Theta_deg"}
        normalized_cols = [col_alias.get(str(col), col) for col in columns]

        with open(filepath, "w") as f:
            f.write("\t".join(str(col) for col in normalized_cols) + "\n")
            for item in items:
                values = table.item(item)["values"]
                keep_text = table.item(item, "text")
                if str(keep_text) not in ("☑", "✓", "True", "1"):
                    continue
                f.write("\t".join(str(val) for val in values) + "\n")

        messagebox.showinfo("Export successful", f"Results exported to:\n{filepath}")
    except Exception as e:
        messagebox.showerror("Export error", f"Error exporting to TXT: {e}")


def on_tree_click(event):
    """Toggle the Keep checkbox when the tree column is clicked."""
    if table is None:
        return
    row_id = table.identify_row(event.y)
    col_id = table.identify_column(event.x)
    if not row_id or col_id != "#0":
        return
    try:
        current_text = table.item(row_id, "text")
        new_text = "☐" if str(current_text) == "☑" else "☑"
        table.item(row_id, text=new_text)
    except Exception:
        pass


def get_latest_chimerax_run_dir():
    """Return the newest timestamped ChimeraX run folder, or the current directory as a fallback."""
    global latest_chimerax_run_dir
    cwd = os.getcwd()
    if latest_chimerax_run_dir and os.path.isdir(latest_chimerax_run_dir):
        return latest_chimerax_run_dir

    candidates = []
    for entry in os.listdir(cwd):
        path = os.path.join(cwd, entry)
        if entry.startswith("ChimeraX_") and os.path.isdir(path):
            candidates.append(path)

    if candidates:
        latest_chimerax_run_dir = max(candidates, key=os.path.getmtime)
        return latest_chimerax_run_dir

    return cwd


def create_chimerax_run_dir():
    """Create and return a timestamped folder for one ChimeraX run."""
    global latest_chimerax_run_dir
    run_dir = os.path.join(os.getcwd(), f"ChimeraX_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    os.makedirs(run_dir, exist_ok=True)
    latest_chimerax_run_dir = run_dir
    return run_dir


def run_chimerax():
    """Run ChimeraX commands per kept table row."""
    try:
        if table is None:
            messagebox.showwarning("No data", "Please compute results first.")
            return
        clear_png_gallery()
        rows = []
        for item in table.get_children():
            if str(table.item(item, "text")) not in ("☑", "✓", "True", "1"):
                continue
            rows.append(table.item(item)["values"])
        if not rows:
            messagebox.showwarning("No selection", "Please keep at least one row to send to ChimeraX.")
            return

        pdb_path = var_pdb_path.get().strip()
        if not pdb_path:
            messagebox.showwarning("Missing PDB", "Please select a PDB model.")
            return

        bin_path = var_chimerax_bin.get().strip() or CHIMERAX_BIN
        if os.path.isdir(bin_path):
            bin_candidate = os.path.join(bin_path, "ChimeraX")
            bin_path = bin_candidate if os.path.exists(bin_candidate) else bin_path
        if not os.path.exists(bin_path):
            messagebox.showerror("ChimeraX not found", f"ChimeraX binary not found at:\n{bin_path}")
            return

        save_chimerax_bin(bin_path)

        try:
            p_size = float(var_p_size.get())
        except Exception:
            p_size = 4.0
        try:
            lp_filter = float(var_lp_filter.get())
        except Exception:
            lp_filter = 15.0
        try:
            n_su = int(float(var_n_su.get()))
        except Exception:
            n_su = 10 if var_mode.get() == "monomers" else 5

        auto_proj_targets = []
        auto_generated_pngs = []
        run_dir = create_chimerax_run_dir()

        for idx, vals in enumerate(rows):
            try:
                mt_type = str(vals[0])
                N_val = float(vals[1])
                S_val = float(vals[2])
                theta_deg = float(vals[3])
                r_su_A = float(vals[4])
                phi_su_deg = float(vals[5])
                delta_c_a = float(vals[6])
                r_pf_A = float(vals[7])
                phi_pf_deg = float(vals[8])
            except Exception:
                continue

            pdb_out_path = os.path.join(run_dir, f"{mt_type}.pdb")
            mrc_out_path = os.path.join(run_dir, f"{mt_type}.mrc")
            png_out_path = os.path.join(run_dir, f"{mt_type}_proj.png")

            cmds = [
                f"open \"{pdb_path}\"",
                f"open \"{pdb_path}\"",
                f"turn 0,1,0 {theta_deg} models #2 coordinateSystem #1",
                f"sym #2 h,{r_su_A},{phi_su_deg},{n_su},{-1*n_su/2} coord #1 center 0,{delta_c_a},0 copies true",
                f"sym #3 h,{r_pf_A},{phi_pf_deg},{N_val} coord #1 center 0,{delta_c_a},0 copies true",
                f"move 0,1,0 {-1*delta_c_a} coord #4 models #4",
                f"save \"{pdb_out_path}\" model #4",
                f"volume new gridSpacing {p_size}",
                f"volume cover #5 atomBox #4",
                f"molmap #4 {lp_filter} onGrid #6",
                f"save \"{mrc_out_path}\" #7",
                "close #1-3",
                "close #5-6"
            ]

            cmd_str = "; ".join(cmds)
            args = [bin_path, "--nogui", "--cmd", cmd_str]
            try:
                for old_path in (mrc_out_path, png_out_path):
                    try:
                        if os.path.exists(old_path):
                            os.remove(old_path)
                    except Exception:
                        pass

                subprocess.Popen(args)
            except Exception as e:
                messagebox.showerror("ChimeraX error", f"Failed to launch ChimeraX for {mt_type}: {e}")
                return

            auto_proj_targets.append(mrc_out_path)

        info_msg = "ChimeraX launched headless (--nogui)."
        proj_success = 0
        proj_errors = []
        display_ok = True
        display_errors = []
        if auto_proj_targets:
            wait_deadline = time.time() + 300.0
            pending = set(auto_proj_targets)
            while pending and time.time() < wait_deadline:
                ready = [p for p in pending if os.path.exists(p) and os.path.getsize(p) > 0]
                for p in ready:
                    pending.remove(p)
                if pending:
                    time.sleep(0.5)
            if pending:
                proj_errors.extend([f"{os.path.basename(p)} not found or empty after wait" for p in sorted(pending)])
            existing = [p for p in auto_proj_targets if os.path.exists(p) and os.path.getsize(p) > 0]
            if existing:
                proj_success, proj_errors_extra, auto_generated_pngs = project_mrc_files(existing, show_message=False)
                if proj_errors_extra:
                    proj_errors.extend(proj_errors_extra)
                if auto_generated_pngs:
                    display_ok, display_errors = refresh_png_gallery(auto_generated_pngs)

        if proj_success:
            info_msg += f"\nAuto-projected {proj_success} .mrc file(s) to PNG."
        if proj_errors:
            info_msg += "\nProjection issues:\n" + "\n".join(str(e) for e in proj_errors)
        if auto_generated_pngs and not display_ok and display_errors:
            info_msg += "\nDisplay issues:\n" + "\n".join(display_errors)

        try:
            refresh_generated_models()
        except Exception:
            pass

        messagebox.showinfo("ChimeraX", info_msg)
    except Exception as e:
        messagebox.showerror("ChimeraX error", str(e))


def project_mrc_files(paths, show_message=True):
    """Project provided .mrc paths along Y, rotate 90°, save as .png."""
    if not paths:
        return 0, [], []

    if not MRCFILE_AVAILABLE or not NP_AVAILABLE:
        missing = []
        if not MRCFILE_AVAILABLE:
            missing.append("mrcfile")
        if not NP_AVAILABLE:
            missing.append("numpy")
        if show_message:
            messagebox.showerror(
                "Missing dependencies",
                "\n".join(
                    ["These packages are required: " + ", ".join(missing), "Install with:", "pip install numpy mrcfile Pillow"]
                ),
            )
        return 0, ["Missing dependencies: " + ", ".join(missing)], []

    save_with_pillow = PIL_AVAILABLE
    save_with_matplotlib = (not save_with_pillow) and MATPLOTLIB_AVAILABLE
    if not save_with_pillow and not save_with_matplotlib:
        if show_message:
            messagebox.showerror(
                "Missing saver",
                "Pillow or matplotlib is required to save PNG files. Install with:\npip install Pillow"
            )
        return 0, ["Missing Pillow/matplotlib for saving PNGs"], []

    success = 0
    errors = []
    generated_pngs = []
    for path in paths:
        try:
            with mrcfile.open(path, permissive=True) as mrc:
                data = mrc.data
            if data is None or getattr(data, "ndim", 0) < 3:
                raise ValueError("MRC data is not 3D")

            proj = np.mean(data, axis=1)  # project along Y
            proj_rot = np.rot90(proj)  # rotate 90°
            proj_rot = np.flipud(proj_rot)  # reflect across X axis

            vmin = np.nanmin(proj_rot)
            vmax = np.nanmax(proj_rot)
            if not np.isfinite(vmin) or not np.isfinite(vmax):
                raise ValueError("Projection contains non-finite values")
            if vmax - vmin < 1e-12:
                norm = np.zeros_like(proj_rot, dtype=np.uint8)
            else:
                norm = ((proj_rot - vmin) / (vmax - vmin) * 255).astype(np.uint8)

            png_path = os.path.splitext(path)[0] + "_proj.png"
            if save_with_pillow:
                Image.fromarray(norm).save(png_path)
            else:
                import matplotlib.pyplot as plt  # local import; safe due to earlier check

                plt.imsave(png_path, norm, cmap="gray")

            success += 1
            generated_pngs.append(png_path)
        except Exception as exc:
            errors.append(f"{os.path.basename(path)}: {exc}")

    if show_message:
        if success:
            message = f"Saved {success} projection(s) to PNG."
            if errors:
                message += "\nErrors:\n" + "\n".join(errors)
            messagebox.showinfo("Projection complete", message)
        elif errors:
            messagebox.showerror("Projection failed", "\n".join(errors))
    return success, errors, generated_pngs


def project_mrc_to_png():
    """Manual entry point: select .mrc maps, project, and save as .png."""
    file_paths = filedialog.askopenfilenames(title="Select .mrc maps", filetypes=[("MRC files", "*.mrc"), ("All files", "*.*")])
    if not file_paths:
        return
    success, errors, pngs = project_mrc_files(list(file_paths), show_message=True)
    if pngs:
        refresh_png_gallery(pngs)


def refresh_png_gallery(png_paths):
    """Display PNG projections in the ChimeraX tab."""
    if gallery_inner_frame is None:
        return False, ["Gallery frame is not initialized"]
    if not PIL_AVAILABLE or ImageTk is None:
        messagebox.showerror(
            "Missing Pillow",
            "Pillow is required to display PNG projections. Install with:\npip install Pillow",
        )
        return False, ["Missing Pillow (ImageTk)"]

    for widget in gallery_inner_frame.winfo_children():
        widget.destroy()
    png_thumbnails.clear()

    if not png_paths:
        return True, []

    count = len(png_paths)
    gallery_inner_frame.update_idletasks()
    available_width = gallery_inner_frame.winfo_width()
    if available_width <= 1:
        try:
            available_width = gallery_inner_frame.master.winfo_width()
        except Exception:
            available_width = 0
    if available_width <= 1:
        available_width = 1200

    min_thumb = 120
    max_thumb = 280
    ideal_cols = max(1, int(available_width / (min_thumb + 12)))
    cols = max(1, min(count, ideal_cols))
    thumb_width = max(min_thumb, min(max_thumb, int(available_width / cols) - 16))
    max_width = thumb_width
    max_height = 260

    errors = []
    for idx, path in enumerate(png_paths):
        try:
            img = Image.open(path)
            img.thumbnail((max_width, max_height))
            photo = ImageTk.PhotoImage(img)
            png_thumbnails.append(photo)

            frame = tk.Frame(gallery_inner_frame, bd=1, relief=tk.SOLID, padx=4, pady=4)
            tk.Label(frame, image=photo).pack()
            tk.Label(frame, text=os.path.basename(path), font=(None, 9)).pack(pady=(4, 0))

            r = idx // cols
            c = idx % cols
            frame.grid(row=r, column=c, padx=6, pady=6, sticky="n")
        except Exception as exc:
            errors.append(f"{os.path.basename(path)}: {exc}")

    if errors:
        for err in errors:
            tk.Label(gallery_inner_frame, text=f"Error loading {err}").grid(sticky="w")
    return len(errors) == 0, errors


def refresh_generated_models():
    """Populate combobox with generated N_S.pdb files from the newest ChimeraX run folder."""
    global generated_combo
    if generated_combo is None:
        return
    try:
        search_dir = get_latest_chimerax_run_dir()
        pdbs = []
        for root, _, files in os.walk(search_dir):
            for name in files:
                if name.lower().endswith(".pdb") and "_" in name:
                    rel_path = os.path.relpath(os.path.join(root, name), os.getcwd())
                    pdbs.append(rel_path)
        pdbs = sorted({p for p in pdbs})
        generated_combo["values"] = pdbs
        if pdbs:
            if var_selected_generated_pdb.get() not in pdbs:
                var_selected_generated_pdb.set(pdbs[0])
        else:
            var_selected_generated_pdb.set("")
    except Exception as exc:
        messagebox.showerror("Refresh error", f"Could not list generated models: {exc}")


def open_selected_generated():
    """Open chosen N_S.pdb (and matching .mrc if present) with viewer defaults."""
    pdb_name = var_selected_generated_pdb.get().strip()
    if not pdb_name:
        messagebox.showwarning("No selection", "Select a generated N_S.pdb from the list or refresh it.")
        return

    pdb_path = pdb_name if os.path.isabs(pdb_name) else os.path.join(os.getcwd(), pdb_name)
    if not os.path.exists(pdb_path):
        messagebox.showerror("Missing file", f"PDB not found: {pdb_path}")
        return

    mrc_path = os.path.splitext(pdb_path)[0] + ".mrc"

    bin_path = var_chimerax_bin.get().strip() or CHIMERAX_BIN
    if os.path.isdir(bin_path):
        bin_candidate = os.path.join(bin_path, "ChimeraX")
        bin_path = bin_candidate if os.path.exists(bin_candidate) else bin_path
    if not os.path.exists(bin_path):
        messagebox.showerror("ChimeraX not found", f"ChimeraX binary not found at:\n{bin_path}")
        return

    save_chimerax_bin(bin_path)

    cmds = [f"open \"{pdb_path}\""]
    mrc_opened = False
    if os.path.exists(mrc_path):
        cmds.append(f"open \"{mrc_path}\"")
        mrc_opened = True
    cmds += [
        "camera ortho",
        "set bgColor dark gray",
        "view all",
        "color bychain",
        "transparency 50",
        "volume #2 level 0.07",
        "color #2 #b2ffff80 models",
    ]

    args = [bin_path, "--cmd", "; ".join(cmds)]
    try:
        subprocess.Popen(args)
        msg = f"Opened {os.path.basename(pdb_path)} in ChimeraX with default view settings."
        if mrc_opened:
            msg += f" Also opened {os.path.basename(mrc_path)}."
        else:
            msg += " (No matching .mrc file found.)"
        messagebox.showinfo("ChimeraX", msg)
    except Exception as exc:
        messagebox.showerror("ChimeraX error", f"Failed to open selection: {exc}")


def clear_png_gallery():
    """Clear PNG projection thumbnails."""
    if gallery_inner_frame is None:
        return
    for widget in gallery_inner_frame.winfo_children():
        widget.destroy()
    png_thumbnails.clear()


def on_mode_change():
    """Update default n_su when switching between monomers and dimers."""
    try:
        if var_mode.get() == "monomers":
            var_n_su.set("10")
        else:
            var_n_su.set("5")
    except Exception:
        pass


def sync_chimerax_mode_from_dimer(*_args):
    """If dimer mode is selected in the Parameters tab, match it in the ChimeraX tab."""
    try:
        if var_use_dimer.get():
            var_mode.set("dimers")
            on_mode_change()
    except Exception:
        pass


def build_gui(root):
    root.title("TubulePy")
    root.geometry("1000x820+30+30")
    root.resizable(True, True)
    try:
        root.update()
        root.lift()
        root.attributes("-topmost", True)
        root.after(500, lambda: root.attributes("-topmost", False))
    except Exception:
        pass

    formula_frm = tk.Frame(root, bg="lightgrey", padx=12, pady=8)
    formula_frm.pack(fill=tk.X)
    tk.Label(
        formula_frm,
        text="Computes helical parameters to generate microtubule N_S references using a PDB model",
        font=("Arial", 12, "italic"),
        bg="lightgrey",
    ).pack()

    notebook = ttk.Notebook(root)
    notebook.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
    tab_main = tk.Frame(notebook)
    tab_chimera = tk.Frame(notebook)
    notebook.add(tab_main, text="Parameters")
    notebook.add(tab_chimera, text="ChimeraX")

    left_frm = tk.LabelFrame(tab_main, text="Parameters", padx=10, pady=10)
    left_frm.pack(side=tk.LEFT, fill=tk.BOTH, padx=(0, 12))

    tk.Label(left_frm, text="N (protofilaments):").grid(row=0, column=0, sticky="w", pady=2)
    tk.Entry(left_frm, textvariable=var_N, width=4).grid(row=0, column=1, sticky="w")

    tk.Label(left_frm, text="S (helical starts):").grid(row=1, column=0, sticky="w", pady=2)
    tk.Entry(left_frm, textvariable=var_S, width=4).grid(row=1, column=1, sticky="w")

    tk.Label(left_frm, text="a (subunit repeat, Å):").grid(row=2, column=0, sticky="w", pady=2)
    tk.Entry(left_frm, textvariable=var_a, width=6).grid(row=2, column=1, sticky="w")

    tk.Label(left_frm, text="r (subunit rise, Å):").grid(row=3, column=0, sticky="w", pady=2)
    r_entry = tk.Entry(left_frm, textvariable=var_r, width=6)
    r_entry.grid(row=3, column=1, sticky="w")
    r_entry.bind("<FocusOut>", on_r_focus_out)

    tk.Label(left_frm, text="δx (separation, Å):").grid(row=4, column=0, sticky="w", pady=2)
    tk.Entry(left_frm, textvariable=var_dx, width=6).grid(row=4, column=1, sticky="w")

    btn_single = tk.Button(left_frm, text="Calculate Single", command=compute_theta_single, width=14)
    btn_single.grid(row=5, column=0, columnspan=2, pady=(8, 4))

    global lbl_result
    lbl_result = tk.Label(left_frm, text="θ = —", font=(None, 10, "bold"), wraplength=150)
    lbl_result.grid(row=6, column=0, columnspan=2, sticky="w", pady=(0, 4))

    global var_use_dimer
    var_use_dimer = tk.BooleanVar(value=False)
    var_use_dimer.trace_add("write", sync_chimerax_mode_from_dimer)
    tk.Checkbutton(left_frm, text="Dimer (double a)", variable=var_use_dimer, bg="LightGrey").grid(
        row=7, column=0, columnspan=2, sticky="w", pady=(0, 8)
    )

    global var_add_refs
    var_add_refs = tk.BooleanVar(value=False)
    tk.Checkbutton(left_frm, text="Add references", variable=var_add_refs, bg="LightGrey").grid(
        row=8, column=0, columnspan=2, sticky="w", pady=(0, 8)
    )

    tk.Frame(left_frm, height=2, bg="gray").grid(row=9, column=0, columnspan=2, sticky="ew", pady=8)

    tk.Label(left_frm, text="N range (min–max):", font=("Arial", 10, "bold")).grid(
        row=10, column=0, columnspan=2, sticky="w", pady=(4, 2)
    )
    tk.Label(left_frm, text="Min:").grid(row=11, column=0, sticky="w", padx=(10, 0))
    tk.Entry(left_frm, textvariable=var_N_min, width=4).grid(row=11, column=1, sticky="w")
    tk.Label(left_frm, text="Max:").grid(row=12, column=0, sticky="w", padx=(10, 0))
    tk.Entry(left_frm, textvariable=var_N_max, width=4).grid(row=12, column=1, sticky="w")

    tk.Label(left_frm, text="S range (min–max):", font=("Arial", 10, "bold")).grid(
        row=13, column=0, columnspan=2, sticky="w", pady=(4, 2)
    )
    tk.Label(left_frm, text="Min:").grid(row=14, column=0, sticky="w", padx=(10, 0))
    tk.Entry(left_frm, textvariable=var_S_min, width=4).grid(row=14, column=1, sticky="w")
    tk.Label(left_frm, text="Max:").grid(row=15, column=0, sticky="w", padx=(10, 0))
    tk.Entry(left_frm, textvariable=var_S_max, width=4).grid(row=15, column=1, sticky="w")

    tk.Checkbutton(left_frm, text="θ filter enabled", variable=var_use_theta_filter, bg="LightGrey").grid(
        row=16, column=0, columnspan=2, sticky="w", pady=(4, 2)
    )
    tk.Label(left_frm, text="θ range (min–max):", font=("Arial", 10, "bold")).grid(
        row=17, column=0, columnspan=2, sticky="w", pady=(2, 2)
    )
    tk.Label(left_frm, text="Min:").grid(row=18, column=0, sticky="w", padx=(10, 0))
    tk.Entry(left_frm, textvariable=var_theta_min, width=4).grid(row=18, column=1, sticky="w")
    tk.Label(left_frm, text="Max:").grid(row=19, column=0, sticky="w", padx=(10, 0))
    tk.Entry(left_frm, textvariable=var_theta_max, width=4).grid(row=19, column=1, sticky="w")

    btn_range = tk.Button(left_frm, text="Calculate Range", command=compute_theta_range, width=14, bg="#e0e0e0")
    btn_range.grid(row=21, column=0, columnspan=2, pady=(4, 4))

    btn_plot = tk.Button(left_frm, text="Plot N ~ θ", command=plot_n_vs_theta, width=14, bg="#d0e8d0")
    btn_plot.grid(row=22, column=0, columnspan=2, pady=(0, 4))

    btn_save = tk.Button(left_frm, text="Save results", command=export_to_txt, width=14, bg="#ffffcc")
    btn_save.grid(row=23, column=0, columnspan=2, pady=(4, 4))

    right_frm = tk.LabelFrame(tab_main, text="Results", padx=10, pady=10)
    right_frm.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    global table
    columns = ("N_S", "N", "S", "θ (°)", "r_su_A", "Phi_su_deg", "Delta_C_A", "r_pf_A", "Phi_pf_deg")
    table = ttk.Treeview(right_frm, columns=columns, height=15, show="tree headings")

    table.column("#0", width=40, anchor=tk.CENTER, stretch=False)
    table.heading("#0", text="Keep")

    table.column("N_S", width=10, anchor=tk.CENTER)
    table.column("N", width=10, anchor=tk.CENTER)
    table.column("S", width=10, anchor=tk.CENTER)
    table.column("θ (°)", width=20, anchor=tk.CENTER)
    table.column("r_su_A", width=20, anchor=tk.CENTER)
    table.column("Phi_su_deg", width=40, anchor=tk.CENTER)
    table.column("Delta_C_A", width=40, anchor=tk.CENTER)
    table.column("r_pf_A", width=40, anchor=tk.CENTER)
    table.column("Phi_pf_deg", width=40, anchor=tk.CENTER)

    table.heading("N_S", text="N_S")
    table.heading("N", text="N")
    table.heading("S", text="S")
    table.heading("θ (°)", text="θ (°)")
    table.heading("r_su_A", text="Z_su (Å)")
    table.heading("Phi_su_deg", text="φ_su (°)")
    table.heading("Delta_C_A", text="δ_radius (Å)")
    table.heading("r_pf_A", text="Z_pf (Å)")
    table.heading("Phi_pf_deg", text="φ_pf (°)")

    scrollbar = ttk.Scrollbar(right_frm, orient=tk.VERTICAL, command=table.yview)
    table.configure(yscroll=scrollbar.set)
    table.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    table.bind("<Button-1>", on_tree_click)

    chx_frm = tk.Frame(tab_chimera, padx=12, pady=12)
    chx_frm.pack(fill=tk.BOTH, expand=True)

    def browse_pdb():
        path = filedialog.askopenfilename(
            title="Select PDB model", filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )
        if path:
            var_pdb_path.set(path)

    def browse_chimerax_bin():
        path = filedialog.askdirectory(
            title="Select ChimeraX directory", initialdir=os.path.dirname(var_chimerax_bin.get() or CHIMERAX_BIN)
        )
        if path:
            candidate = os.path.join(path, "ChimeraX")
            resolved = candidate if os.path.exists(candidate) else path
            var_chimerax_bin.set(resolved)
            save_chimerax_bin(resolved)

    row_idx = 0
    tk.Label(chx_frm, text="ChimeraX binary:").grid(row=row_idx, column=0, sticky="w", pady=2)
    tk.Entry(chx_frm, textvariable=var_chimerax_bin, width=50).grid(
        row=row_idx, column=1, columnspan=2, sticky="we", pady=2
    )
    tk.Button(chx_frm, text="Browse", command=browse_chimerax_bin).grid(row=row_idx, column=3, padx=4)
    row_idx += 1

    tk.Label(chx_frm, text="PDB model:").grid(row=row_idx, column=0, sticky="w", pady=2)
    tk.Entry(chx_frm, textvariable=var_pdb_path, width=50).grid(row=row_idx, column=1, columnspan=2, sticky="we", pady=2)
    tk.Button(chx_frm, text="Browse", command=browse_pdb).grid(row=row_idx, column=3, padx=4)
    row_idx += 1

    tk.Label(chx_frm, text="Mode:").grid(row=row_idx, column=1, sticky="e", pady=2)
    tk.Radiobutton(chx_frm, text="Monomers", variable=var_mode, value="monomers", command=on_mode_change).grid(
        row=row_idx, column=2, sticky="w"
    )
    tk.Radiobutton(chx_frm, text="Dimers", variable=var_mode, value="dimers", command=on_mode_change).grid(
        row=row_idx, column=3, sticky="w"
    )
    row_idx += 1

    controls_frm = tk.Frame(chx_frm)
    controls_frm.grid(row=row_idx, column=0, columnspan=7, sticky="ew", pady=4)

    tk.Label(controls_frm, text="Molecules (n_su):").grid(row=0, column=0, sticky="e", padx=4)
    tk.Entry(controls_frm, textvariable=var_n_su, width=6).grid(row=0, column=1, sticky="w", padx=4)
    tk.Label(controls_frm, text="Pixel size (Å):").grid(row=0, column=2, sticky="e", padx=4)
    tk.Entry(controls_frm, textvariable=var_p_size, width=6).grid(row=0, column=3, sticky="w", padx=4)
    tk.Label(controls_frm, text="Low-pass filter (Å):").grid(row=0, column=4, sticky="e", padx=4)
    tk.Entry(controls_frm, textvariable=var_lp_filter, width=6).grid(row=0, column=5, sticky="w", padx=4)
    tk.Button(controls_frm, text="Run ChimeraX", command=run_chimerax, bg="#d0e8d0", width=14).grid(
        row=0, column=6, sticky="e", padx=6
    )
    for i in range(7):
        controls_frm.columnconfigure(i, weight=1)
    row_idx += 1

    # Picker for generated N_S models in current working directory
    tk.Label(chx_frm, text="Generated N_S model:").grid(row=row_idx, column=0, sticky="w", pady=2)
    global generated_combo
    generated_combo = ttk.Combobox(chx_frm, textvariable=var_selected_generated_pdb, width=40, state="readonly")
    generated_combo.grid(row=row_idx, column=1, columnspan=2, sticky="we", pady=2)
    tk.Button(chx_frm, text="Refresh", command=refresh_generated_models).grid(row=row_idx, column=3, padx=4)
    tk.Button(chx_frm, text="Open in ChimeraX", command=open_selected_generated, bg="#d0e8d0").grid(
        row=row_idx, column=4, padx=4
    )
    row_idx += 1

    gallery_container = tk.LabelFrame(chx_frm, text="PNG projections", padx=8, pady=8)
    gallery_container.grid(row=row_idx, column=0, columnspan=7, sticky="nsew")

    canvas = tk.Canvas(gallery_container, height=260)
    vscrollbar = ttk.Scrollbar(gallery_container, orient=tk.VERTICAL, command=canvas.yview)
    hscrollbar = ttk.Scrollbar(gallery_container, orient=tk.HORIZONTAL, command=canvas.xview)
    canvas.configure(yscrollcommand=vscrollbar.set, xscrollcommand=hscrollbar.set)
    vscrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    hscrollbar.pack(side=tk.BOTTOM, fill=tk.X)
    canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    global gallery_inner_frame
    gallery_inner_frame = tk.Frame(canvas)
    canvas.create_window((0, 0), window=gallery_inner_frame, anchor="nw")

    def _on_frame_configure(event):
        canvas.configure(scrollregion=canvas.bbox("all"))

    gallery_inner_frame.bind("<Configure>", _on_frame_configure)
    gallery_container.columnconfigure(0, weight=1)

    chx_frm.rowconfigure(row_idx, weight=1)

    for i in range(7):
        chx_frm.columnconfigure(i, weight=1)

    refresh_generated_models()


if __name__ == "__main__":
    root = tk.Tk()

    def _auto_messagebox(kind, title, message, timeout_ms=5000):
        """Show a simple dialog that auto-closes after timeout."""
        win = tk.Toplevel(root)
        win.title(title)
        win.transient(root)
        win.grab_set()
        bg = "white"
        fg = {"info": "black", "warning": "darkorange", "error": "darkred"}.get(kind, "black")
        header = tk.Label(win, text=title, font=(None, 11, "bold"), fg=fg, bg=bg, anchor="w", justify="left")
        header.pack(fill="x", padx=12, pady=(10, 4))
        body = tk.Message(win, text=message, width=480, bg=bg, anchor="w", justify="left")
        body.pack(fill="both", padx=12, pady=4)
        tk.Button(win, text="OK", command=win.destroy).pack(pady=(0, 10))
        win.configure(bg=bg)
        try:
            win.after(timeout_ms, win.destroy)
        except Exception:
            pass
        try:
            win.lift()
            win.focus_force()
        except Exception:
            pass

    def _wrap_msg(kind):
        def _inner(title, message, **kwargs):
            return _auto_messagebox(kind, title, message, timeout_ms=kwargs.get("timeout", 5000))
        return _inner

    messagebox.showinfo = _wrap_msg("info")
    messagebox.showwarning = _wrap_msg("warning")
    messagebox.showerror = _wrap_msg("error")

    var_N = tk.StringVar(value="13")
    var_S = tk.StringVar(value="3")
    var_a = tk.StringVar(value="40.60")
    var_r_exact = 9.36923077
    var_r = tk.StringVar(value=f"{var_r_exact:.2f}")
    var_dx = tk.StringVar(value="51.40")

    var_N_min = tk.StringVar(value="10")
    var_N_max = tk.StringVar(value="16")
    var_S_min = tk.DoubleVar(value=2.0)
    var_S_max = tk.DoubleVar(value=4.0)

    var_use_theta_filter = tk.BooleanVar(value=True)
    var_theta_min = tk.StringVar(value="-3.5")
    var_theta_max = tk.StringVar(value="3.5")

    lbl_result = None
    table = None
    var_add_refs = None

    var_map_path = tk.StringVar(value="")
    var_pdb_path = tk.StringVar(value=os.path.join(os.getcwd(), "PDB_models", "Alpha.pdb"))
    var_mode = tk.StringVar(value="monomers")
    var_n_su = tk.StringVar(value="10")
    var_p_size = tk.StringVar(value="4")
    var_lp_filter = tk.StringVar(value="15")
    var_chimerax_bin = tk.StringVar(value=CHIMERAX_BIN)
    var_selected_generated_pdb = tk.StringVar(value="")
    gallery_inner_frame = None

    build_gui(root)

    import sys

    print(f"GUI started with Python: {sys.executable}")
    print(f"Matplotlib available: {MATPLOTLIB_AVAILABLE}")

    root.mainloop()

"""Written by Denis Chrétien, using VS Code IDE and Tkinter library.
   Designed with AI assistance.
   Last updated: March 2026.
   Contact: denis.chretien@univ-rennes.fr"""
