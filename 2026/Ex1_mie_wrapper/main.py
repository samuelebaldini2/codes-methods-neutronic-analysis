import miescat

# reference data
m_real = 1.33
m_img = 0.01
radius = 1.0
wavelength = 0.55

print("==================================================================")
print("compute_mie_scattering")
print("==================================================================")

result = miescat.compute_mie_scattering(
    m_real,
    m_img,
    radius,
    wavelength,
)

for key, value in result.items():
    print(f"{key:<35} : {value}")



