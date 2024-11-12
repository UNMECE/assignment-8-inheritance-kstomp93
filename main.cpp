#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

// Base class Field
class Field {
protected:
    double *value; 

public:
    // Default constructor 
    Field() {
        value = new double[3]{0.0, 0.0, 0.0};
    }

    // Parameterized constructor
    Field(double x, double y, double z) {
        value = new double[3]{x, y, z};
    }

    // Destructor
    virtual ~Field() {
        delete[] value;
    }

    // Copy constructor 
    Field(const Field &f) {
        value = new double[3];
        for (int i = 0; i < 3; ++i)
            value[i] = f.value[i];
    }

    // Print magnitude function
    virtual void printMagnitude() const {
        std::cout << "Components: (" << value[0] << ", " << value[1] << ", " << value[2] << ")" << std::endl;
    }

};

// Derived class Electric_Field
class Electric_Field : public Field {
private:
    double calculated_E; 

public:
    // Default constructor
    Electric_Field() : Field() {}

    // Parameterized constructor 
    Electric_Field(double ex, double ey, double ez) : Field(ex, ey, ez) {}

    // Getters
    double getX() const { return value[0]; }
    double getY() const { return value[1]; }
    double getZ() const { return value[2]; }

    // Setters
    void setX(double ex) { value[0] = ex; }
    void setY(double ey) { value[1] = ey; }
    void setZ(double ez) { value[2] = ez; }

    // Calculate magnitude 
    double calculateMagnitude() const {
        return sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
    }

    // Calculate electric 
    void calculateElectricField(double Q, double r) {
        const double epsilon0 = 8.854e-12; // Permittivity of free space
        calculated_E = Q / (4 * M_PI * epsilon0 * r * r);
    }

    // Overload the '+' operator 
    Electric_Field operator+(const Electric_Field &other) const {
        return Electric_Field(value[0] + other.value[0], value[1] + other.value[1], value[2] + other.value[2]);
    }

    // Print calculated electric field
    void printCalculatedField() const {
        std::cout << "Calculated Electric Field: " << calculated_E << " N/C" << std::endl;
    }

    // Overloaded '<<' operator
    friend std::ostream &operator<<(std::ostream &out, const Electric_Field &e) {
        out << "Electric Field Components: (" << e.value[0] << ", " << e.value[1] << ", " << e.value[2] << ")";
        return out;
    }
};

// Derived class Magnetic_Field
class Magnetic_Field : public Field {
private:
    double calculated_B; 

public:
    // Default constructor
    Magnetic_Field() : Field() {}

    // Parameterized constructor 
    Magnetic_Field(double bx, double by, double bz) : Field(bx, by, bz) {}

    // Getters
    double getX() const { return value[0]; }
    double getY() const { return value[1]; }
    double getZ() const { return value[2]; }

    // Setters
    void setX(double bx) { value[0] = bx; }
    void setY(double by) { value[1] = by; }
    void setZ(double bz) { value[2] = bz; }

    // Calculate magnitude 
    double calculateMagnitude() const {
        return sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
    }

    // Calculate magnetic field using Ampere's Law
    void calculateMagneticField(double I, double r) {
        const double mu0 = 4 * M_PI * 1e-7; // Permeability of free space
        calculated_B = (mu0 * I) / (2 * M_PI * r);
    }

    // Overload the '+' operator 
    Magnetic_Field operator+(const Magnetic_Field &other) const {
        return Magnetic_Field(value[0] + other.value[0], value[1] + other.value[1], value[2] + other.value[2]);
    }

    // Print calculated magnetic field
    void printCalculatedField() const {
        std::cout << "Calculated Magnetic Field: " << calculated_B << " T" << std::endl;
    }

    // Overloaded '<<' operator
    friend std::ostream &operator<<(std::ostream &out, const Magnetic_Field &m) {
        out << "Magnetic Field Components: (" << m.value[0] << ", " << m.value[1] << ", " << m.value[2] << ")";
        return out;
    }
};

int main() {
    // Create electric and magnetic field objects with specific components
    Electric_Field E1(1.0, 1e5, 1e3);
    Magnetic_Field B1(0.2, 0.5, 0.8);

    // Modify components using setters
    E1.setX(3.2);
    E1.setY(4.5);
    E1.setZ(6.7);
    B1.setX(0.3);
    B1.setY(0.4);
    B1.setZ(0.5);

    // Print field components using printMagnitude
    E1.printMagnitude();
    B1.printMagnitude();

    // Calculate fields using Gauss's and Ampere's laws
    E1.calculateElectricField(1e-6, 0.1); 
    B1.calculateMagneticField(10, 0.05); 
    
    // Print calculated field values
    E1.printCalculatedField();
    B1.printCalculatedField();

    // Demonstrate overloaded '+' operator
    Electric_Field E2(2.0, 2e5, 2e3);
    Electric_Field E_sum = E1 + E2;
    std::cout << "Sum of Electric Fields: " << E_sum << std::endl;

    Magnetic_Field B2(0.3, 0.6, 0.9);
    Magnetic_Field B_sum = B1 + B2;
    std::cout << "Sum of Magnetic Fields: " << B_sum << std::endl;

    // Demonstrate overloaded '<<' operator
    std::cout << "E1: " << E1 << std::endl;
    std::cout << "B1: " << B1 << std::endl;

    return 0;
}
