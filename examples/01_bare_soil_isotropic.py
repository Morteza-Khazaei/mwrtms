"""Example: Isotropic bare soil backscatter."""

from mwrtms import mwRTMs


def main() -> None:
    result = mwRTMs.compute_soil_backscatter(
        model="aiem",
        frequency_ghz=5.4,
        incident_angle_deg=40.0,
        rms_height_cm=1.0,
        correlation_length_cm=5.0,
        correlation_type="exponential",
        soil_moisture=0.25,
        clay_fraction=0.3,
        sand_fraction=0.5,
    )

    print("AIEM Backscatter Results:")
    print(f"  σ⁰_HH = {result['hh']:.2f} dB")
    print(f"  σ⁰_VV = {result['vv']:.2f} dB")
    print(f"  σ⁰_HV = {result['hv']:.2f} dB")


if __name__ == "__main__":
    main()
