"""Example: Anisotropic bare soil (tilled rows)."""

from mwrtms import mwRTMs


def main() -> None:
    result = mwRTMs.compute_soil_backscatter(
        model="aiem",
        frequency_ghz=5.4,
        incident_angle_deg=40.0,
        rms_height_cm=1.5,
        correlation_length_cm=15.0,
        correlation_length_y_cm=3.0,
        correlation_type="exponential",
        anisotropic=True,
        soil_moisture=0.20,
    )

    print("Anisotropic (tilled soil):")
    print(f"  σ⁰_VV = {result['vv']:.2f} dB")
    print(f"  Anisotropy ratio = {15.0/3.0:.1f}")


if __name__ == "__main__":
    main()
