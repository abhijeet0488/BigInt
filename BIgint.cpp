// Corrected and Optimized C++ program for BigInt operations
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <cmath> // For log10 in count_digits

using namespace std;

class BigInt {
private:
    // Store number in a large base to reduce operations.
    static const int base = 1000000000;
    static const int base_digits = 9;
    vector<int> digits;
    bool is_negative = false;

    // --- Private Helper Enum for internal constructor ---
    enum class PrivateConstructor { E };
    BigInt(PrivateConstructor); // Private constructor for internal use

    // --- Helper Functions ---
    void trim(); // Removes leading zeros
    static void divide_by_2(BigInt& a); // Fast division by 2
    static BigInt karatsuba_multiply(const BigInt& a, const BigInt& b);

public:
    // --- Constructors ---
    BigInt(long long n = 0);
    BigInt(const string& s);
    
    // --- Copy Semantics ---
    BigInt(const BigInt& other);
    BigInt& operator=(const BigInt& other);

    // --- Move Semantics ---
    BigInt(BigInt&& other) noexcept;
    BigInt& operator=(BigInt&& other) noexcept;

    // --- Operator Overloading ---
    BigInt& operator=(long long n);
    BigInt operator-() const;
    friend bool operator==(const BigInt& a, const BigInt& b);
    friend bool operator<(const BigInt& a, const BigInt& b);
    friend BigInt operator+(BigInt a, const BigInt& b);
    friend BigInt operator-(BigInt a, const BigInt& b);
    friend BigInt operator*(const BigInt& a, const BigInt& b);
    friend BigInt operator/(const BigInt& a, const BigInt& b);
    friend BigInt operator%(const BigInt& a, const BigInt& b);
    BigInt& operator+=(const BigInt& other);
    BigInt& operator-=(const BigInt& other);
    BigInt& operator*=(const BigInt& other);
    BigInt& operator/=(const BigInt& other);
    BigInt& operator%=(const BigInt& other);
    BigInt& operator++(); BigInt operator++(int);
    BigInt& operator--(); BigInt operator--(int);

    // --- I/O Operations ---
    friend istream& operator>>(istream& in, BigInt& a);
    friend ostream& operator<<(ostream& out, const BigInt& a);

    // --- Other Functions ---
    bool isZero() const;
    friend BigInt abs(const BigInt& a);
    friend size_t count_digits(const BigInt& a); // ADDED

    // --- Mathematical Functions ---
    friend BigInt gcd(BigInt a, BigInt b);
    friend BigInt power(BigInt base, BigInt exp);
    friend BigInt modular_power(BigInt base, BigInt exp, const BigInt& mod);
    friend BigInt Factorial(int n);
    friend BigInt sqrt(const BigInt& a); // ADDED
    friend BigInt NthFibonacci(int n); // ADDED
    friend BigInt NthCatalan(int n); // ADDED
};

// --- Private Constructor for internal use ---
BigInt::BigInt(PrivateConstructor) {}

// --- Helper Implementation ---
void BigInt::trim() {
    while (digits.size() > 1 && digits.back() == 0) {
        digits.pop_back();
    }
    if (digits.empty() || (digits.size() == 1 && digits[0] == 0)) {
        is_negative = false;
    }
}

void BigInt::divide_by_2(BigInt& a) {
    int carry = 0;
    for (int i = a.digits.size() - 1; i >= 0; i--) {
        long long cur = a.digits[i] + (long long)carry * base;
        a.digits[i] = cur / 2;
        carry = cur % 2;
    }
    a.trim();
}

// --- Constructor Implementation ---
BigInt::BigInt(long long n) {
    if (n == 0) { digits.push_back(0); return; }
    if (n < 0) { is_negative = true; n = -n; }
    while (n > 0) {
        digits.push_back(n % base);
        n /= base;
    }
}

BigInt::BigInt(const string& s) {
    if (s.empty() || s == "0" || s == "-0") { digits.push_back(0); return; }
    int start = 0;
    if (s[0] == '-') { is_negative = true; start = 1; }
    else if (s[0] == '+') { start = 1; }
    for (int i = s.length(); i > start; i -= base_digits) {
        long long chunk = 0;
        int len = max(start, i - base_digits);
        for (int j = len; j < i; ++j) {
            if (!isdigit(s[j])) throw invalid_argument("Invalid character in number string");
            chunk = chunk * 10 + (s[j] - '0');
        }
        digits.push_back(chunk);
    }
    trim();
}

// --- Copy Semantics Implementation ---
BigInt::BigInt(const BigInt& other) : digits(other.digits), is_negative(other.is_negative) {}

BigInt& BigInt::operator=(const BigInt& other) {
    if (this != &other) {
        digits = other.digits;
        is_negative = other.is_negative;
    }
    return *this;
}

// --- Move Semantics Implementation ---
BigInt::BigInt(BigInt&& other) noexcept : digits(move(other.digits)), is_negative(other.is_negative) {}

BigInt& BigInt::operator=(BigInt&& other) noexcept {
    if (this != &other) {
        digits = move(other.digits);
        is_negative = other.is_negative;
    }
    return *this;
}

// --- Operator Implementation ---
BigInt& BigInt::operator=(long long n) {
    digits.clear();
    is_negative = false;
    if (n == 0) { digits.push_back(0); }
    else {
        if (n < 0) { is_negative = true; n = -n; }
        while (n > 0) { digits.push_back(n % base); n /= base; }
    }
    return *this;
}

BigInt BigInt::operator-() const {
    BigInt result = *this;
    if (!result.isZero()) result.is_negative = !result.is_negative;
    return result;
}

bool operator<(const BigInt& a, const BigInt& b) {
    if (a.is_negative != b.is_negative) return a.is_negative;
    if (a.is_negative) return abs(b) < abs(a);
    if (a.digits.size() != b.digits.size()) return a.digits.size() < b.digits.size();
    for (int i = a.digits.size() - 1; i >= 0; --i) {
        if (a.digits[i] != b.digits[i]) return a.digits[i] < b.digits[i];
    }
    return false;
}
bool operator==(const BigInt& a, const BigInt& b) { return a.digits == b.digits && a.is_negative == b.is_negative; }
bool operator!=(const BigInt& a, const BigInt& b) { return !(a == b); }
bool operator>(const BigInt& a, const BigInt& b) { return b < a; }
bool operator<=(const BigInt& a, const BigInt& b) { return !(b < a); }
bool operator>=(const BigInt& a, const BigInt& b) { return !(a < b); }

BigInt& BigInt::operator+=(const BigInt& other) {
    if (is_negative == other.is_negative) {
        int carry = 0;
        for (size_t i = 0; i < max(digits.size(), other.digits.size()) || carry; ++i) {
            if (i == digits.size()) digits.push_back(0);
            long long current = digits[i] + carry + (i < other.digits.size() ? other.digits[i] : 0);
            digits[i] = current % base;
            carry = current / base;
        }
    } else *this -= (-other);
    trim();
    return *this;
}
BigInt operator+(BigInt a, const BigInt& b) { return a += b; }

BigInt& BigInt::operator-=(const BigInt& other) {
     if (is_negative == other.is_negative) {
        if (abs(*this) >= abs(other)) {
            int borrow = 0;
            for (size_t i = 0; i < other.digits.size() || borrow; ++i) {
                long long current = (long long)digits[i] - borrow - (i < other.digits.size() ? other.digits[i] : 0);
                digits[i] = current < 0 ? current + base : current;
                borrow = current < 0 ? 1 : 0;
            }
        } else *this = -(other - *this);
    } else *this += (-other);
    trim();
    return *this;
}
BigInt operator-(BigInt a, const BigInt& b) { return a -= b; }

BigInt operator*(const BigInt& a, const BigInt& b) {
    if (min(a.digits.size(), b.digits.size()) < 32) { // Threshold for Karatsuba
        BigInt result;
        if (a.isZero() || b.isZero()) return result;
        result.digits.resize(a.digits.size() + b.digits.size());
        for (size_t i = 0; i < a.digits.size(); ++i) {
            for (size_t j = 0, carry = 0; j < b.digits.size() || carry > 0; ++j) {
                long long current = result.digits[i + j] + a.digits[i] * (j < b.digits.size() ? (long long)b.digits[j] : 0) + carry;
                result.digits[i + j] = current % BigInt::base;
                carry = current / BigInt::base;
            }
        }
        result.is_negative = a.is_negative != b.is_negative;
        result.trim();
        return result;
    }
    return BigInt::karatsuba_multiply(a, b);
}

BigInt BigInt::karatsuba_multiply(const BigInt& a, const BigInt& b) {
    size_t n = max(a.digits.size(), b.digits.size());
    if (n < 32) return a * b;
    n = (n + 1) / 2;
    BigInt a_lo(PrivateConstructor::E), a_hi(PrivateConstructor::E);
    BigInt b_lo(PrivateConstructor::E), b_hi(PrivateConstructor::E);
    a_lo.digits.assign(a.digits.begin(), a.digits.begin() + min(n, a.digits.size()));
    if (a.digits.size() > n) a_hi.digits.assign(a.digits.begin() + n, a.digits.end());
    b_lo.digits.assign(b.digits.begin(), b.digits.begin() + min(n, b.digits.size()));
    if (b.digits.size() > n) b_hi.digits.assign(b.digits.begin() + n, b.digits.end());
    a_lo.trim(); a_hi.trim(); b_lo.trim(); b_hi.trim();
    BigInt p1 = karatsuba_multiply(a_hi, b_hi);
    BigInt p2 = karatsuba_multiply(a_lo, b_lo);
    BigInt p3 = karatsuba_multiply(a_hi + a_lo, b_hi + b_lo);
    BigInt p3_sub = p3 - p1 - p2;
    p1.digits.insert(p1.digits.begin(), 2 * n, 0);
    p3_sub.digits.insert(p3_sub.digits.begin(), n, 0);
    BigInt result = p1 + p2 + p3_sub;
    result.is_negative = a.is_negative != b.is_negative;
    result.trim();
    return result;
}

BigInt operator/(const BigInt& a, const BigInt& b) {
    if (b.isZero()) throw runtime_error("Arithmetic Error: Division by zero");
    BigInt abs_a = abs(a), abs_b = abs(b);
    if (abs_a < abs_b) return BigInt(0);
    BigInt low(1), high = abs_a, res(0);
    while (low <= high) {
        BigInt mid_diff = high - low;
        BigInt::divide_by_2(mid_diff);
        BigInt mid = low + mid_diff;
        if (mid.isZero()) break;
        if (abs_b * mid <= abs_a) {
            res = mid;
            low = mid + 1;
        } else high = mid - 1;
    }
    res.is_negative = a.is_negative != b.is_negative;
    res.trim();
    return res;
}

BigInt& BigInt::operator*=(const BigInt& other) { *this = *this * other; return *this; }
BigInt& BigInt::operator/=(const BigInt& other) { *this = *this / other; return *this; }
BigInt operator%(const BigInt& a, const BigInt& b) { return a - (a / b) * b; }
BigInt& BigInt::operator%=(const BigInt& other) { *this = *this % other; return *this; }

BigInt& BigInt::operator++() { *this += 1; return *this; }
BigInt BigInt::operator++(int) { BigInt temp = *this; ++(*this); return temp; }
BigInt& BigInt::operator--() { *this -= 1; return *this; }
BigInt BigInt::operator--(int) { BigInt temp = *this; --(*this); return temp; }

istream& operator>>(istream& in, BigInt& a) { string s; in >> s; a = BigInt(s); return in; }
ostream& operator<<(ostream& out, const BigInt& a) {
    if (a.isZero()) return out << 0;
    if (a.is_negative) out << "-";
    out << a.digits.back();
    for (int i = a.digits.size() - 2; i >= 0; --i)
        out << setfill('0') << setw(BigInt::base_digits) << a.digits[i];
    return out;
}

bool BigInt::isZero() const { return digits.size() == 1 && digits[0] == 0; }
BigInt abs(const BigInt& a) { BigInt result = a; result.is_negative = false; return result; }

// --- ADDED: count_digits function ---
size_t count_digits(const BigInt& a) {
    if (a.isZero()) return 1;
    size_t count = (a.digits.size() - 1) * BigInt::base_digits;
    int last_digit = a.digits.back();
    if (last_digit > 0) {
        count += floor(log10(last_digit)) + 1;
    }
    return count;
}


// --- Mathematical Function Implementations ---
BigInt Factorial(int n) { if (n < 0) throw invalid_argument("Factorial of negative number is not defined."); BigInt f(1); for (int i = 2; i <= n; i++) f *= i; return f; }
BigInt gcd(BigInt a, BigInt b) { a = abs(a); b = abs(b); while (!b.isZero()) { a %= b; swap(a, b); } return a; }

BigInt power(BigInt base, BigInt exp) {
    if (exp < 0) throw invalid_argument("Exponent must be non-negative.");
    BigInt res(1);
    while (!exp.isZero()) {
        if (exp.digits[0] & 1) res *= base;
        base *= base;
        BigInt::divide_by_2(exp);
    }
    return res;
}

BigInt modular_power(BigInt base, BigInt exp, const BigInt& mod) {
    if (exp < 0) throw invalid_argument("Exponent must be non-negative.");
    BigInt res(1);
    base %= mod;
    while (!exp.isZero()) {
        if (exp.digits[0] & 1) res = (res * base) % mod;
        base = (base * base) % mod;
        BigInt::divide_by_2(exp);
    }
    return res;
}

// --- ADDED: sqrt function ---
BigInt sqrt(const BigInt& a) {
    if (a.is_negative) throw invalid_argument("sqrt of negative number is not defined.");
    if (a.isZero()) return BigInt(0);
    BigInt low(1), high = a, res(0);
    while(low <= high) {
        BigInt mid = low + (high - low) / 2;
        if(mid.isZero()) { low = 1; continue; }
        BigInt mid_sq = mid * mid;
        if (mid_sq == a) return mid;
        if (mid_sq < a) {
            res = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return res;
}

// --- ADDED: NthFibonacci function ---
BigInt NthFibonacci(int n) {
    if (n < 0) throw invalid_argument("Fibonacci of negative number is not defined.");
    if (n == 0) return BigInt(0);
    if (n <= 2) return BigInt(1);
    BigInt a(1), b(1);
    for (int i = 3; i <= n; ++i) {
        BigInt c = a + b;
        a = b;
        b = c;
    }
    return b;
}

// --- ADDED: NthCatalan function ---
BigInt NthCatalan(int n) {
    if (n < 0) throw invalid_argument("Catalan of negative number is not defined.");
    // Formula: C_n = (1/(n+1)) * (2n choose n) = (2n)! / ((n+1)! * n!)
    BigInt C = Factorial(2 * n) / (Factorial(n + 1) * Factorial(n));
    return C;
}


// Driver code
int main() {
    cout << "--- Basic Tests ---\n";
    BigInt first("12345678901234567890");
    cout << "First Number: " << first << endl;
    cout << "Number of digits in first: " << count_digits(first) << endl;
    
    cout << "\n--- Square Root Test ---\n";
    BigInt to_sqrt("100000000000000000000");
    cout << "Square root of " << to_sqrt << " is " << sqrt(to_sqrt) << endl;

    cout << "\n--- Fibonacci Test ---\n";
    cout << "The 100th Fibonacci number is: " << NthFibonacci(100) << endl;

    cout << "\n--- Catalan Test ---\n";
    cout << "The 50th Catalan number is: " << NthCatalan(50) << endl;
    
    cout << "\n--- Factorial Test ---\n";
    cout << "Factorial of 50 = " << Factorial(50) << endl;

    return 0;
}