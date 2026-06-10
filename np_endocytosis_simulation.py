import numpy as np
import matplotlib.pyplot as plt

# 参数
Rw = 48.6
Dw = 2 * Rw
D = np.linspace(1, 200, 5000)
X = D / Dw

# ===================== CIE 计算 =====================
mbo_cie = -0.1
a_cie = Rw * mbo_cie
y_cie = (X**2 - (1 + a_cie * X)**2) / X**3
y_cie[y_cie < 0] = 0
# 归一化到 1
y_cie_norm = y_cie / np.max(y_cie)
# 找峰值
idx_max_cie = np.argmax(y_cie_norm)
D_max_cie = D[idx_max_cie]
y_max_cie = y_cie_norm[idx_max_cie]

# ===================== CME 计算 =====================
mbo_cme = -0.025
a_cme = Rw * mbo_cme
y_cme = (X**2 - (1 + a_cme * X)**2) / X**3
y_cme[y_cme < 0] = 0
# 归一化到 1
y_cme_norm = y_cme / np.max(y_cme)
# 找峰值
idx_max_cme = np.argmax(y_cme_norm)
D_max_cme = D[idx_max_cme]
y_max_cme = y_cme_norm[idx_max_cme]

# ===================== 绘图：CIE =====================
plt.figure(figsize=(8, 5))
color_cie = '#5797D0'

plt.plot(D, y_cie_norm, color=color_cie, linewidth=2.5)
# 峰值垂线
plt.plot([D_max_cie, D_max_cie], [0, y_max_cie], color=color_cie, linestyle='--', linewidth=1.5)
# 峰值标注
plt.text(D_max_cie + 14, y_max_cie*0.5, f'{D_max_cie:.1f} nm',
         color=color_cie, ha='center', fontsize=11, weight='bold')

# 0 刻度线
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)

# 坐标轴设置
plt.xlabel('NP Diameter (nm)', fontsize=12)
plt.ylabel('Endocytosis rate', fontsize=12)
plt.xlim(0, 200)
plt.ylim(0, 1.05)  # 归一化后上限 1
plt.yticks([])     # 不显示 y 轴数字
plt.grid(False)    # 关闭网格
plt.tight_layout()

# ===================== 绘图：CME =====================
plt.figure(figsize=(8, 5))
color_cme = '#BF699D'

plt.plot(D, y_cme_norm, color=color_cme, linewidth=2.5)
# 峰值垂线
plt.plot([D_max_cme, D_max_cme], [0, y_max_cme], color=color_cme, linestyle='--', linewidth=1.5)
# 峰值标注
plt.text(D_max_cme + 12, y_max_cme*0.5, f'{D_max_cme:.1f} nm',
         color=color_cme, ha='center', fontsize=11, weight='bold')

# 0 刻度线
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)

# 坐标轴设置
plt.xlabel('NP Diameter (nm)', fontsize=12)
plt.ylabel('Endocytosis rate', fontsize=12)
plt.xlim(0, 200)
plt.ylim(0, 1.05)
plt.yticks([])
plt.grid(False)
plt.tight_layout()

# 显示所有图
plt.show()

# 输出数值
print(f"CIE 峰值直径 = {D_max_cie:.1f} nm")
print(f"CME 峰值直径 = {D_max_cme:.1f} nm")
