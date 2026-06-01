import numpy as np

def test_dispersion():
    npts = 100
    kx = np.linspace(-0.5, 0.5, npts)
    ky = np.linspace(-0.5, 0.5, npts)
    kz = np.linspace(-0.5, 0.5, npts)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    # BCC
    E_bcc = -8.0 * np.cos(2*np.pi*KX) * np.cos(2*np.pi*KY) * np.cos(2*np.pi*KZ)
    print(f"BCC: {E_bcc.min():.4f} to {E_bcc.max():.4f}")

    # FCC
    E_fcc = -4.0 * (np.cos(2*np.pi*KX)*np.cos(2*np.pi*KY) + np.cos(2*np.pi*KY)*np.cos(2*np.pi*KZ) + np.cos(2*np.pi*KZ)*np.cos(2*np.pi*KX))
    print(f"FCC: {E_fcc.min():.4f} to {E_fcc.max():.4f}")

    # Hexagonal
    E_hex = -2.0 * (np.cos(4*np.pi*KX) + 2.0 * np.cos(2*np.pi*KX) * np.cos(2*np.pi * np.sqrt(3) * KY)) - np.cos(2*np.pi*KZ)
    print(f"Hexagonal: {E_hex.min():.4f} to {E_hex.max():.4f}")

if __name__ == "__main__":
    test_dispersion()
