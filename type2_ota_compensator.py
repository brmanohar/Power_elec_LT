"""
Type II Compensator Design using OTA
Based on TI Application Report SLVA662 and LTspice design methodology
(as shown in the "Simulating a Type 2 Compensator" slide)

Design equations from slide:
    boost  = PM - PS - 90
    G      = 10^(-Gfc/20)
    k      = tan((boost/2 + 45) * pi/180)      [Type II uses boost/2, not boost/4]
    fp     = fc * k
    fz     = fc / k
    a      = sqrt((fc^2/fp^2) + 1)
    b      = sqrt((fz^2/fc^2) + 1)
    R2     = (a/b) * (fp*G) * (Rlower+Rupper) / ((fp-fz)*Rlower*gm)
    C1     = 1 / (2*pi*R2*fz)
    C2     = (Rlower*gm / (2*pi*fp*G*(Rlower+Rupper))) * (b/a)

References:
    - TI SLVA662: "Demystifying Type II and Type III Compensators Using OpAmp and OTA"
    - LTspice .VAR equations from design slide
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# =============================================================================
# DESIGN INPUTS
# =============================================================================

# --- Extract from plant (measured or simulated) ---
Gfc_dB  = -20.0     # Plant gain at crossover frequency [dB]
PS_deg  = -70.0     # Plant phase lag at crossover [degrees]

# --- Design targets ---
fc      = 1e3       # Target crossover frequency [Hz]
PM      = 70.0      # Desired phase margin [degrees]

# --- System / hardware parameters ---
Vout    = 12.0      # Output voltage [V]
Ibias   = 250e-6    # Bias current [A]
Vref    = 2.5       # Reference voltage [V]
gm      = 10e-6     # OTA transconductance [S]

# =============================================================================
# DERIVED PARAMETERS
# =============================================================================

Rlower = Vref / Ibias
Rupper = (Vout - Vref) / Ibias

print("=" * 60)
print("  Type II OTA Compensator Design")
print("=" * 60)
print(f"\n[Inputs]")
print(f"  Plant gain at fc   : {Gfc_dB} dB")
print(f"  Plant phase at fc  : {PS_deg}°")
print(f"  Crossover freq fc  : {fc/1e3:.1f} kHz")
print(f"  Phase margin PM    : {PM}°")
print(f"  Vout               : {Vout} V")
print(f"  Vref               : {Vref} V")
print(f"  Ibias              : {Ibias*1e6:.0f} µA")
print(f"  gm (OTA)           : {gm*1e6:.0f} µS")

print(f"\n[Resistor Divider]")
print(f"  Rlower = Vref/Ibias             = {Rlower/1e3:.2f} kΩ")
print(f"  Rupper = (Vout-Vref)/Ibias      = {Rupper/1e3:.2f} kΩ")

# =============================================================================
# STEP 1: Required compensator gain at crossover
# =============================================================================

G_lin = 10 ** (-Gfc_dB / 20)
print(f"\n[Step 1] Required compensator gain at fc")
print(f"  G = 10^(-Gfc/20) = {G_lin:.4f}  ({-Gfc_dB:.1f} dB)")

# =============================================================================
# STEP 2: Required phase boost
# Type II: boost = PM - PS - 90  (same formula, but Type II max boost = 90°)
# =============================================================================

boost_deg = PM - PS_deg - 90.0
print(f"\n[Step 2] Required phase boost")
print(f"  boost = PM - PS - 90 = {PM} - ({PS_deg}) - 90 = {boost_deg:.1f}°")

if boost_deg >= 90:
    print(f"  ⚠ WARNING: boost={boost_deg:.1f}° ≥ 90°. Type II max is ~90°.")
    print(f"    Consider using Type III compensator instead.")

# =============================================================================
# STEP 3: Pole/zero placement
# Type II uses a SINGLE pole/zero pair (asymmetric k formula: boost/2 + 45)
# =============================================================================

k  = np.tan((boost_deg / 2 + 45) * np.pi / 180)   # Type II: boost/2
fp = fc * k
fz = fc / k

print(f"\n[Step 3] Pole/Zero placement")
print(f"  k  = tan((boost/2 + 45)·π/180) = {k:.6f}")
print(f"  fp = fc · k  = {fp:.2f} Hz")
print(f"  fz = fc / k  = {fz:.2f} Hz")

# =============================================================================
# STEP 4: Correction factors
# =============================================================================

a_val = np.sqrt((fc / fp) ** 2 + 1)
b_val = np.sqrt((fz / fc) ** 2 + 1)

print(f"\n[Step 4] Magnitude correction factors")
print(f"  a = sqrt((fc/fp)²+1) = {a_val:.6f}")
print(f"  b = sqrt((fz/fc)²+1) = {b_val:.6f}")

# =============================================================================
# STEP 5: Component calculation  (from slide .VAR equations)
#   R2 = (a/b) * fp*G*(Rlower+Rupper) / ((fp-fz)*Rlower*gm)
#   C1 = 1 / (2π·R2·fz)
#   C2 = (Rlower·gm / (2π·fp·G·(Rlower+Rupper))) · (b/a)
# =============================================================================

R2 = (a_val / b_val) * (fp * G_lin * (Rlower + Rupper)) / ((fp - fz) * Rlower * gm)
C1 = 1.0 / (2 * np.pi * R2 * fz)
C2 = (Rlower * gm / (2 * np.pi * fp * G_lin * (Rlower + Rupper))) * (b_val / a_val)

print(f"\n[Step 5] Passive Component Values")
print(f"  R2 = {R2:.2f} Ω   ({R2/1e3:.4f} kΩ)")
print(f"  C1 = {C1*1e12:.2f} pF  ({C1*1e9:.4f} nF)")
print(f"  C2 = {C2*1e12:.2f} pF  ({C2*1e9:.4f} nF)")

# =============================================================================
# VERIFICATION: Frequency response
# Type II OTA transfer function (SLVA662 eq. 12/13):
#   H(s) = -[R4·gm/(R1+R4)] · (1 + sR2·C1) / [s(C1+C2)(1 + sR2·C1‖C2)]
# where R1=Rupper, R4=Rlower, C3→C2 (only one pole cap in Type II)
# =============================================================================

freqs = np.logspace(1, 5, 3000)
s     = 1j * 2 * np.pi * freqs

C1C2_ser = (C1 * C2) / (C1 + C2)   # C1 || C2 (series combination for HF pole)

num = 1 + s * R2 * C1
den = s * (C1 + C2) * (1 + s * R2 * C1C2_ser)

H = (Rlower * gm / (Rupper + Rlower)) * (num / den)

gain_dB     = 20 * np.log10(np.abs(H))
phase_raw   = np.angle(H) * 180 / np.pi
phase_boost = phase_raw + 90.0      # shift: integrator baseline (−90°) → 0°

idx_fc      = np.argmin(np.abs(freqs - fc))
gain_at_fc  = gain_dB[idx_fc]
boost_at_fc = phase_boost[idx_fc]

print(f"\n[Verification at fc = {fc/1e3:.1f} kHz]")
print(f"  Compensator gain   : {gain_at_fc:.2f} dB  (target: {-Gfc_dB:.1f} dB)")
print(f"  Phase boost        : {boost_at_fc:.1f}°   (target: {boost_deg:.1f}°)")

# =============================================================================
# PLOT
# =============================================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
fig.suptitle(
    f"Type II OTA Compensator – Bode Plot\n"
    f"fc = {fc/1e3:.0f} kHz,  PM = {PM:.0f}°,  Boost = {boost_deg:.0f}°",
    fontsize=13
)

# --- Gain ---
ax1.semilogx(freqs, gain_dB, 'k', linewidth=2, label='|H(f)|')
ax1.axvline(fc, color='navy',  linestyle='--', linewidth=1.2, label=f'fc = {fc/1e3:.1f} kHz')
ax1.axvline(fz, color='green', linestyle=':',  linewidth=1.2, label=f'fz = {fz:.0f} Hz')
ax1.axvline(fp, color='red',   linestyle=':',  linewidth=1.2, label=f'fp = {fp:.0f} Hz')
ax1.axhline(gain_at_fc, color='navy', linestyle=':', linewidth=0.8, alpha=0.5)
ax1.annotate(
    f'{gain_at_fc:.1f} dB @ fc',
    xy=(fc, gain_at_fc),
    xytext=(fc * 2.5, gain_at_fc + 4),
    color='navy', fontsize=9,
    arrowprops=dict(arrowstyle='->', color='navy')
)
ax1.set_ylabel("Gain / dB", fontsize=11)
ax1.legend(fontsize=8, loc='lower left')
ax1.grid(True, which='both', alpha=0.3)

# --- Phase Boost ---
ax2.semilogx(freqs, phase_boost, 'r', linewidth=2, label='Phase boost')
ax2.axvline(fc, color='navy',  linestyle='--', linewidth=1.2)
ax2.axvline(fz, color='green', linestyle=':',  linewidth=1.2)
ax2.axvline(fp, color='red',   linestyle=':',  linewidth=1.2)
ax2.axhline(0,          color='gray',  linestyle='-',  linewidth=0.7, alpha=0.5)
ax2.axhline(boost_at_fc, color='navy', linestyle=':', linewidth=0.8, alpha=0.5)
ax2.annotate(
    f'{boost_at_fc:.0f}° boost @ fc',
    xy=(fc, boost_at_fc),
    xytext=(fc * 2.5, boost_at_fc - 8),
    color='navy', fontsize=9,
    arrowprops=dict(arrowstyle='->', color='navy')
)
ax2.set_ylabel("Phase boost / degrees\n(above integrator −90°)", fontsize=10)
ax2.set_xlabel("Frequency / Hz", fontsize=11)
ax2.grid(True, which='both', alpha=0.3)
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(
    lambda x, _: f'{x/1e3:.0f}k' if x >= 1000 else f'{x:.0f}'
))

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/type2_ota_bode.png", dpi=150, bbox_inches='tight')
plt.show()
print("\n[Plot saved] type2_ota_bode.png")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 60)
print("  FINAL COMPONENT SUMMARY")
print("=" * 60)
print(f"  {'Component':<12} {'Value':>16}")
print(f"  {'-'*30}")
print(f"  {'Rupper':<12} {Rupper/1e3:>13.2f} kΩ")
print(f"  {'Rlower':<12} {Rlower/1e3:>13.2f} kΩ")
print(f"  {'R2':<12} {R2/1e3:>13.4f} kΩ")
print(f"  {'C1':<12} {C1*1e12:>13.2f} pF")
print(f"  {'C2':<12} {C2*1e12:>13.2f} pF")
print(f"  {'fz':<12} {fz:>13.2f} Hz")
print(f"  {'fp':<12} {fp:>13.2f} Hz")
print(f"  {'Boost':<12} {boost_deg:>13.1f}°")
print("=" * 60)
