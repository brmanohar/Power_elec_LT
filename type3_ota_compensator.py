"""
Type III Compensator Design using OTA
Based on TI Application Report SLVA662 and LTspice design methodology
(as shown in the slide with moderate phase boost approach)

Design approach: symmetric pole/zero placement around crossover frequency
    fz1 = fz2 = fc / sqrt(Kf)
    fp1 = fp2 = fc * sqrt(Kf)
    where Kf = (tan(boost/4 + 45) * pi/180)^2

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
Gfc_dB   = -20.0    # Plant gain at crossover frequency [dB]
PS_deg   = -135.0   # Plant phase lag at crossover frequency [degrees]

# --- Design targets ---
fc       = 1e3      # Target crossover frequency [Hz]
PM       = 50.0     # Desired phase margin [degrees]

# --- System / hardware parameters ---
Vout     = 19.0     # Output voltage [V]
Ibias    = 250e-6   # Bias current [A]
Vref     = 2.5      # Reference voltage [V]
gm       = 10e-6    # OTA transconductance [S]

# =============================================================================
# DERIVED PARAMETERS
# =============================================================================

Rlower = Vref / Ibias                    # Lower feedback resistor [Ω]
Rupper = (Vout - Vref) / Ibias           # Upper feedback resistor [Ω]

print("=" * 60)
print("  Type III OTA Compensator Design")
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

G_lin = 10 ** (-Gfc_dB / 20)   # Compensator must provide gain = 1/|plant gain|
print(f"\n[Step 1] Required compensator gain at fc")
print(f"  G = 10^(-Gfc/20) = 10^({-Gfc_dB}/20) = {G_lin:.4f}  ({-Gfc_dB:.1f} dB)")

# =============================================================================
# STEP 2: Required phase boost
# =============================================================================

boost_deg = PM - PS_deg - 90.0
print(f"\n[Step 2] Required phase boost")
print(f"  boost = PM - PS - 90 = {PM} - ({PS_deg}) - 90 = {boost_deg:.1f}°")

# =============================================================================
# STEP 3: Kf factor and symmetric pole/zero placement
# =============================================================================

Kf = (np.tan((boost_deg / 4 + 45) * np.pi / 180)) ** 2

fz1 = fc / np.sqrt(Kf)
fz2 = fc / np.sqrt(Kf)   # symmetric: fz1 = fz2
fp1 = fc * np.sqrt(Kf)
fp2 = fc * np.sqrt(Kf)   # symmetric: fp1 = fp2

print(f"\n[Step 3] Pole/Zero placement (symmetric)")
print(f"  Kf  = (tan(boost/4 + 45)·π/180)² = {Kf:.6f}")
print(f"  fz1 = fc/sqrt(Kf) = {fz1:.2f} Hz")
print(f"  fz2 = fc/sqrt(Kf) = {fz2:.2f} Hz")
print(f"  fp1 = fc·sqrt(Kf) = {fp1:.2f} Hz")
print(f"  fp2 = fc·sqrt(Kf) = {fp2:.2f} Hz")

# =============================================================================
# STEP 4: Correction factors a, b, c, d (magnitude at crossover)
# =============================================================================
# From TI SLVA662 equations (51)-(55) and the LTspice .VAR definitions

a_val = np.sqrt((fc / fp2) ** 2 + 1)
b_val = np.sqrt((fc / fp1) ** 2 + 1)
c_val = np.sqrt((fz1 / fc) ** 2 + 1)
d_val = np.sqrt((fc / fz2) ** 2 + 1)

aa = (a_val * b_val) / (c_val * d_val)

print(f"\n[Step 4] Magnitude correction factors")
print(f"  a = sqrt((fc/fp2)²+1) = {a_val:.6f}")
print(f"  b = sqrt((fc/fp1)²+1) = {b_val:.6f}")
print(f"  c = sqrt((fz1/fc)²+1) = {c_val:.6f}")
print(f"  d = sqrt((fc/fz2)²+1) = {d_val:.6f}")
print(f"  aa = (a·b)/(c·d)       = {aa:.6f}")

# =============================================================================
# STEP 5: Component calculation
# From SLVA662 equations (56)-(64) combined with LTspice .VAR approach
#
#   bb = G * fp1 * (Rupper + Rlower) / (Rlower * gm * (fp1 - fz1))
#   R2 = aa * bb
#   R3 = (cc) / (dd)
#   C1 = 1 / (2π·fz1·R2)
#   C2 = 1 / (2π·fz2·(R3 + Rupper))   [note: R1=Rupper in OTA topology]
#   C3 = 1 / (2π·fz2·(Rupper + R3))   via C3 = C1/(2π·C1·R2·fp2 - 1)
# =============================================================================

# R2
bb = (G_lin * fp1 * (Rupper + Rlower)) / (Rlower * gm * (fp1 - fz1))
R2 = aa * bb

# R3  (from LTspice .VAR cc/dd)
cc = (Rupper ** 2) * fz2 - Rupper * Rlower * (fp1 - fz2)
dd = (fp1 - fz2) * (Rupper + Rlower)
R3 = cc / dd

# Capacitors
C1 = 1.0 / (2 * np.pi * fz1 * R2)
C2 = 1.0 / (2 * np.pi * fz2 * (Rupper + R3))   # Rupper plays role of R1
C3 = 1.0 / (2 * np.pi * fz2 * (Rupper + R3))   # same as C2 in symmetric case
# Recalculate C3 from SLVA662 eq. 64
C3 = C1 / (2 * np.pi * C1 * R2 * fp2 - 1)

print(f"\n[Step 5] Passive Component Values")
print(f"  bb = G·fp1·(Rupper+Rlower)/(Rlower·gm·(fp1-fz1)) = {bb:.4f}")
print(f"  R2 = aa·bb  = {R2:.2f} Ω   ({R2/1e3:.4f} kΩ)")
print(f"  R3 = cc/dd  = {R3:.2f} Ω   ({R3/1e3:.4f} kΩ)")
print(f"  C1 = 1/(2π·fz1·R2)           = {C1*1e12:.2f} pF  ({C1*1e9:.4f} nF)")
print(f"  C2 = 1/(2π·fz2·(Rupper+R3))  = {C2*1e12:.2f} pF  ({C2*1e9:.4f} nF)")
print(f"  C3 = C1/(2π·C1·R2·fp2-1)     = {C3*1e12:.2f} pF  ({C3*1e9:.4f} nF)")

# =============================================================================
# VERIFICATION: Compute frequency response of the designed compensator
# =============================================================================

freqs = np.logspace(1, 5, 2000)   # 10 Hz to 100 kHz
s = 1j * 2 * np.pi * freqs

# Type III OTA transfer function (SLVA662 eq. 43 / 45)
#   H(s) = -[R4·gm/(R4+R1)] · [1 + sC2(R1+R3)] / [1 + sC2(R4‖R1 + R3)]
#                             · [1 + sR2·C1]     / [s(C1+C3)(1 + sR2·C1‖C3)]
#
# Using Rupper as R1, Rlower as R4

R_par = (Rupper * Rlower) / (Rupper + Rlower)   # Rupper || Rlower

num1 = 1 + s * C2 * (Rupper + R3)
den1 = 1 + s * C2 * (R_par + R3)

num2 = 1 + s * R2 * C1
den2 = s * (C1 + C3) * (1 + s * R2 * C1 * C3 / (C1 + C3))

H = (Rlower * gm / (Rupper + Rlower)) * (num1 / den1) * (num2 / den2)

gain_dB   = 20 * np.log10(np.abs(H))
phase_raw = np.angle(H) * 180 / np.pi          # raw phase: −90° to +90° range
# Phase boost = how much above a pure integrator (−90°).
# A pure integrator sits at −90°; adding +90° maps that to 0°,
# so the boost curve is centred around 0 and peaks positively.
phase_boost = phase_raw + 90.0                  # shift: integrator baseline → 0°

# Find values at crossover
idx_fc = np.argmin(np.abs(freqs - fc))
gain_at_fc  = gain_dB[idx_fc]
boost_at_fc = phase_boost[idx_fc]

print(f"\n[Verification at fc = {fc/1e3:.1f} kHz]")
print(f"  Compensator gain   : {gain_at_fc:.2f} dB  (target: {-Gfc_dB:.1f} dB)")
print(f"  Phase (raw)        : {phase_raw[idx_fc]:.2f}°")
print(f"  Phase boost        : {boost_at_fc:.1f}°  (target: {boost_deg:.1f}°)")

# =============================================================================
# PLOT
# =============================================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
fig.suptitle("Type III OTA Compensator – Bode Plot\n"
             f"fc={fc/1e3:.0f} kHz, PM={PM}°, boost={boost_deg:.0f}°", fontsize=13)

# Gain plot
ax1.semilogx(freqs, gain_dB, 'k', linewidth=2)
ax1.axvline(fc, color='navy', linestyle='--', linewidth=1.2, label=f'fc = {fc/1e3:.1f} kHz')
ax1.axvline(fz1, color='green', linestyle=':', linewidth=1, label=f'fz = {fz1:.0f} Hz')
ax1.axvline(fp1, color='red', linestyle=':', linewidth=1, label=f'fp = {fp1:.0f} Hz')
ax1.axhline(gain_at_fc, color='navy', linestyle=':', linewidth=0.8, alpha=0.5)
ax1.annotate(f'{gain_at_fc:.1f} dB', xy=(fc, gain_at_fc),
             xytext=(fc * 2, gain_at_fc + 3), color='navy', fontsize=9)
ax1.set_ylabel("Gain / dB", fontsize=11)
ax1.legend(fontsize=8, loc='lower left')
ax1.grid(True, which='both', alpha=0.3)
ax1.set_ylim([-20, max(gain_dB) + 10])

# Phase plot  (phase boost relative to integrator baseline)
ax2.semilogx(freqs, phase_boost, 'r', linewidth=2)
ax2.axvline(fc, color='navy', linestyle='--', linewidth=1.2)
ax2.axvline(fz1, color='green', linestyle=':', linewidth=1)
ax2.axvline(fp1, color='red', linestyle=':', linewidth=1)
ax2.axhline(0, color='gray', linestyle='-', linewidth=0.7, alpha=0.5)  # integrator reference
ax2.axhline(boost_at_fc, color='navy', linestyle=':', linewidth=0.8, alpha=0.5)
ax2.annotate(f'{boost_at_fc:.0f}° boost\n@ fc', xy=(fc, boost_at_fc),
             xytext=(fc * 2.5, boost_at_fc - 15), color='navy', fontsize=9,
             arrowprops=dict(arrowstyle='->', color='navy'))
ax2.set_ylabel("Phase boost / degrees\n(above integrator −90°)", fontsize=10)
ax2.set_xlabel("Frequency / Hz", fontsize=11)
ax2.grid(True, which='both', alpha=0.3)
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(
    lambda x, _: f'{x/1e3:.0f}k' if x >= 1000 else f'{x:.0f}'))

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/type3_ota_bode.png", dpi=150, bbox_inches='tight')
plt.show()
print("\n[Plot saved] type3_ota_bode.png")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("\n" + "=" * 60)
print("  FINAL COMPONENT SUMMARY")
print("=" * 60)
print(f"  {'Component':<12} {'Value':>15}")
print(f"  {'-'*30}")
print(f"  {'Rupper':<12} {Rupper/1e3:>12.2f} kΩ")
print(f"  {'Rlower':<12} {Rlower/1e3:>12.2f} kΩ")
print(f"  {'R2':<12} {R2/1e3:>12.4f} kΩ")
print(f"  {'R3':<12} {R3:>12.4f} Ω")
print(f"  {'C1':<12} {C1*1e12:>12.2f} pF")
print(f"  {'C2':<12} {C2*1e12:>12.2f} pF")
print(f"  {'C3':<12} {C3*1e12:>12.2f} pF")
print(f"  {'fz1=fz2':<12} {fz1:>12.2f} Hz")
print(f"  {'fp1=fp2':<12} {fp1:>12.2f} Hz")
print(f"  {'Boost':<12} {boost_deg:>12.1f}°")
print("=" * 60)
