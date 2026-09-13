import miescat

# reference data
m_real = 1.02239
m_img = -0.0119
radius = 2
wavelength = 0.5

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



