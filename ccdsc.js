// by astro / astronomy487
// thank u alex / alexolas

class ccdsc {
	static CCDSC = class {
		constructor(x, i, j, ij, e, ie, je, ije) {
			if (ije != undefined) { //constructor 1: given all parts
				this.real = x;
				this.i = i;
				this.j = j;
				this.ij = ij;
				this.e = e;
				this.ie = ie;
				this.je = je;
				this.ije = ije;
			} else if (typeof x == "object") { //constructor 1.5: given all parts as object
				this.real = this.i = this.j = this.ij = this.e = this.ie = this.je = this.ije = 0;
				for (let k of Object.keys(x))
					this[k] = x[k];
			} else if (typeof x == "number") { //constructor 2: given real x
				this.real = x;
				this.i = this.j = this.ij = this.e = this.ie = this.je = this.ije = 0;
			} else if (typeof x == "string") { //constructor 3: given text
				this.real = this.i = this.j = this.ij = this.e = this.ie = this.je = this.ije = 0;
				let parts = x.replaceAll("ϵ", "e").replaceAll(" ", "").replaceAll("-","+-").split("+").filter(x => x.length);
				for (let component of parts) {
					let num = parseFloat(component);
					if (Number.isNaN(num)) num = component.startsWith("-") ? -1 : 1; //implied 1 or -1
					let iParts = (component.match(/i/g) || []).length;
					while (iParts >= 2) {
						iParts -= 2;
						num *= -1;
					}
					let jParts = (component.match(/j/g) || []).length;
					jParts %= 2;
					let eParts = (component.match(/e/g) || []).length;
					if (eParts >= 2) break; //any term with ee literally doesn't matter
					let property = (iParts?"i":"") + (jParts?"j":"") + (eParts?"e":"");
					property ||= "real";
					this[property] += num;
				}
			} else throw new Error("Couldn't construct a CCDSC");
			for (let property of ["real","i","j","ij","e","ie","je","ije"])
				if (Number.isNaN(this[property]) || Math.abs(this[property]) == Infinity)
					throw new Error("CCDSC number includes infinities or NaNs");
		}
		toString(digits = 5) {
			let txt = "";
			for (let property of ["real","i","j","ij","e","ie","je","ije"]) {
				if (this[property]) {
					let sign = (txt.length&&this[property]>0?"+":"") + (this[property]<0?"-":"");
					let numberText = rounded(Math.abs(this[property]));
					if (numberText == 0) continue;
					if (property != "real" && numberText == 1) numberText = "";
					let unit = property=="real"?"":property.replace("e","ϵ");
					txt += sign + " " + numberText + unit + " ";
				}
			}
			if (!txt.length) txt = "0";
			return txt.trim();
			function rounded(x) {
				if (digits) x = Math.round(x * Math.pow(10, digits))/Math.pow(10, digits);
				return x;
			}
		}
		toTex(digits = 5) {
			return this.toString(digits).replaceAll("ϵ", "\\epsilon ");
		}
		toHtml(digits = 5) {
			return this.toString(digits).replaceAll("i", "<var>i</var>").replaceAll("j", "<var>j</var>").replaceAll("ϵ", "<var>ϵ</var>");
		}
		plus(other) {
			return new ccdsc.CCDSC(
				this.real + other.real,
				this.i + other.i,
				this.j + other.j,
				this.ij + other.ij,
				this.e + other.e,
				this.ie + other.ie,
				this.je + other.je,
				this.ije + other.ije
			);
		}
		minus(other) {
			return new ccdsc.CCDSC(
				this.real - other.real,
				this.i - other.i,
				this.j - other.j,
				this.ij - other.ij,
				this.e - other.e,
				this.ie - other.ie,
				this.je - other.je,
				this.ije - other.ije
			);
		}
		times(other) {
			let [a, b, c, d, e, f, g, h, α, β, γ, δ, ε, ϝ, ζ, η] = [this.real, this.i, this.e, this.j, this.ie, this.ij, this.je, this.ije, other.real, other.i, other.e, other.j, other.ie, other.ij, other.je, other.ije];
			return new ccdsc.CCDSC(
				a*α-b*β+d*δ-f*ϝ,
				a*β+b*α+f*δ+d*ϝ,
				a*δ-b*ϝ+d*α-f*β,
				a*ϝ+b*δ+d*β+f*α,
				a*γ+c*α-b*ε-e*β+d*ζ+g*δ-f*η-h*ϝ,
				a*ε+b*γ+c*β+d*η+e*α+f*ζ+g*ϝ+h*δ,
				a*ζ-b*η+c*δ+d*γ-e*ϝ-f*ε+g*α-h*β,
				a*η+b*ζ+c*ϝ+d*ε+e*δ+f*γ+g*β+h*α
			);
		}
		exp() {
			let answer = new ccdsc.CCDSC(Math.exp(this.real));
			answer = answer.times(new ccdsc.CCDSC({real: Math.cos(this.i), i: Math.sin(this.i)}));
			answer = answer.times(new ccdsc.CCDSC({real: Math.cosh(this.j), j: Math.sinh(this.j)}));
			answer = answer.times(new ccdsc.CCDSC({real: Math.cos(this.ij), ij: Math.sin(this.ij)}));
			answer = answer.times(new ccdsc.CCDSC({real: 1, e: this.e}));
			answer = answer.times(new ccdsc.CCDSC({real: 1, ie: this.ie}));
			answer = answer.times(new ccdsc.CCDSC({real: 1, je: this.je}));
			answer = answer.times(new ccdsc.CCDSC({real: 1, ije: this.ije}));
			return answer;
		}
		scale(scalar) {
			return new ccdsc.CCDSC(
				this.real * scalar,
				this.i * scalar,
				this.j * scalar,
				this.ij * scalar,
				this.e * scalar,
				this.ie * scalar,
				this.je * scalar,
				this.ije * scalar
			);
		}
		magnitude() {
			return Math.sqrt([this.real, this.i, this.j, this.ij].reduce((a,b)=>a+b*b,0))
		}
		over(other) {
			function rref(r){for(var f=r.length,n=r[0].length,e=0,o=0;o<f;o++){if(n<=e)return;for(var t=o;0===r[t][e];)if(t++,f===t&&(t=o,n===++e))return;var a=r[t],i=r[o];r[t]=i,r[o]=a;for(var u=r[o][e],v=0;v<n;v++)r[o][v]/=u;for(t=0;t<f;t++)if(t!==o){u=r[t][e];for(v=0;v<n;v++)r[t][v]-=u*r[o][v]}e++}return r}
			let [A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P] = [
				this.real, this.i, this.j, this.ij, this.e, this.ie, this.je, this.ije,
				other.real, other.i, other.j, other.ij, other.e, other.ie, other.je, other.ije
			];
			let augMatrix = [
				[I, -J, K, -L, 0, 0, 0, 0, A],
				[J, I, L, K, 0, 0, 0, 0, B],
				[K, -L, I, -J, 0, 0, 0, 0, C],
				[L, K, J, I, 0, 0, 0, 0, D],
				[M, -N, O, -P, I, -J, K, -L, E],
				[N, M, P, O, J, I, L, K, F],
				[O, -P, M, -N, K, -L, I, -J, G],
				[P, O, N, M, L, K, J, I, H]
			];
			
			//det(left 8x8) = J**8 + 4*J**6*K**2 - 4*J**6*L**2 + 4*I**2*J**6 - 16*I*J**5*K*L + 6*J**4*K**4 - 4*J**4*K**2*L**2 + 4*I**2*J**4*K**2 + 6*J**4*L**4 - 4*I**2*J**4*L**2 + 6*I**4*J**4 - 32*I*J**3*K**3*L + 32*I*J**3*K*L**3 - 32*I**3*J**3*K*L + 4*J**2*K**6 + 4*J**2*K**4*L**2 - 4*I**2*J**2*K**4 - 4*J**2*K**2*L**4 + 88*I**2*J**2*K**2*L**2 - 4*I**4*J**2*K**2 - 4*J**2*L**6 - 4*I**2*J**2*L**4 + 4*I**4*J**2*L**2 + 4*I**6*J**2 - 16*I*J*K**5*L - 32*I*J*K**3*L**3 + 32*I**3*J*K**3*L - 16*I*J*K*L**5 - 32*I**3*J*K*L**3 - 16*I**5*J*K*L + K**8 + 4*K**6*L**2 - 4*I**2*K**6 + 6*K**4*L**4 - 4*I**2*K**4*L**2 + 6*I**4*K**4 + 4*K**2*L**6 + 4*I**2*K**2*L**4 - 4*I**4*K**2*L**2 - 4*I**6*K**2 + L**8 + 4*I**2*L**6 + 6*I**4*L**4 + 4*I**6*L**2 + I**8;
			
			augMatrix = rref(augMatrix);
			if (!augMatrix) throw new Error("Couldnt divide");
			//if left 8x8 is identity, then read last column for
			let leftIdentity = true;
			for (let i = 0; i < 8; i++)
				for (let j = 0; j < 8; j++)
					if (augMatrix[i][j] != (i==j)?1:0) leftIdentity = false;
			if (!leftIdentity) throw new Error("Couldn't divide");
			return new ccdsc.CCDSC(...augMatrix.map(x => x[8]));
		}
		ln() {
			function complexSplitComplexLn(w, x, y, z) {
				let ogInput = new ccdsc.CCDSC(w,x,y,z,0,0,0,0);
				//https://www.quora.com/Can-logarithms-and-exponentials-be-defined-and-calculated-for-split-complex-numbers-in-a-similar-way-as-with-standard-complex-numbers
				// this = w + xi + yj + zij
				// a = w + xi
				// b = y + zi
				// ln(a+bj) = 1/2 ln (aa - bb) + j/2 ln ((a+b)/(a-b))
				// ln(w+xi+yj+zij) = 1/2 ln (aa-bb) + j/2 ln ((a+b)(a+b)/(a-b)(a+b))
				// ln(w+xi+yj+zij) = 1/2 ln (aa-bb) + j/2 ln ((aa+2ab+bb)/(aa-bb))
				// aa = (w+xi)(w+xi) = ww-xx + 2wxi
				// bb = yy-zz + 2yzi
				let aa = [w*w-x*x, 2*w*x];
				let bb = [y*y-z*z, 2*y*z];
				let aa_minus_bb = [aa[0]-bb[0], aa[1]-bb[1]];
				//aa+2ab+bb = ww-xx + 2wxi + yy-zz + 2yzi + 2(wy-xz+wzi+xiy)
				//aa+2ab+bb = aa + bb + 2(wy-xz+wzi+xiy)
				let aa_2ab_bb = [aa[0]+bb[0]+2*w*y-2*x*z, aa[1]+bb[1]+2*w*z+2*x*y];
				// ln_aa_minus_bb = [real, im]
				// ln_aa_2ab_bb_over_aa_minus_bb = [real, im]
				let ln_aa_minus_bb = complexLn(...aa_minus_bb);
				let ln_aa_2ab_bb_over_aa_minus_bb = complexLn(...complexDiv(...aa_2ab_bb, ...aa_minus_bb));
				let list = [...ln_aa_minus_bb.map(x => x/2), ...ln_aa_2ab_bb_over_aa_minus_bb.map(x => x/2)];
				for (let entry of list) if (Number.isNaN(entry) || Math.abs(entry) == Infinity) throw new Error("CSC logarithm was about to give " + JSON.stringify(list));
				if (isCorrect(list)) return list;
					else list[1] += Math.PI;
				if (isCorrect(list)) return list;
					else list[3] += Math.PI;
				if (isCorrect(list)) return list;
					else list[1] -= Math.PI;
				if (isCorrect(list)) return list;
					else throw new Error("Got wrong branch of logarithm, but adding i*pi didnt fix it. sorry");
				
				function isCorrect(list) {
					let a = new ccdsc.CCDSC(...list, 0, 0, 0, 0).exp();
					let valid = a.minus(ogInput).magnitude() < 10e-10;
					return valid;
				}
				
				function complexDiv(a, b, c, d) {
					//(a+bi)/(c+di) = (a+bi)(c-di)/(c+di)(c-di) = (ac+bd-adi+bic)/(cc+dd)
					let cc_dd = c*c+d*d;
					return [(a*c+b*d)/cc_dd, (b*c-a*d)/cc_dd];
				}
			}
			function complexLn(a, b) {
				return [Math.log(Math.sqrt(a*a+b*b)), Math.atan2(b, a)];
			}
			try {
				if (!(this.j || this.ij || this.e || this.ie || this.je || this.ije)) return new ccdsc.CCDSC(...complexLn(this.real, this.i), 0, 0, 0, 0, 0, 0);
				if (!(this.e || this.ie || this.je || this.ije)) return new ccdsc.CCDSC(...complexSplitComplexLn(this.real, this.i, this.j, this.ij), 0, 0, 0, 0);
				let ln_A = complexSplitComplexLn(this.real, this.i, this.j, this.ij);
				let B_over_A = new CCDSC(this.e, this.ie, this.je, this.ije, 0, 0, 0, 0).over(new CCDSC(this.real, this.i, this.j, this.ij, 0, 0, 0, 0));
				return new ccdsc.CCDSC(
					...ln_A,
					B_over_A.real, B_over_A.i, B_over_A.j, B_over_A.ij
				);
			} catch(e) {
				throw new Error("Logarithm of " + this.toString() + " is undefined");
			}
		}
		equals(other, tolerance = 10e-10) {
			for (let property of ["real","i","j","ij","e","ie","je","ije"])
				if (Math.abs(this[property] - other[property]) > tolerance) return false;
			return true;
		}
	}
	static make(x) {
		try {
			return new this.CCDSC(x);
		} catch(e) {
			throw new Error("Couldn't make a CCDSC number using input "+JSON.stringify(x));
		}
	}
	static ZERO = new this.CCDSC(0);
	static ONE = new this.CCDSC(1);
	static NEGATIVE_ONE = new this.CCDSC(-1);
	static I = new this.CCDSC("i");
	static NEGATIVE_I = new this.CCDSC("-i");
	static NEGATIVE_HALF_I = new this.CCDSC("-0.5i");
	static HALF_I = new this.CCDSC("0.5i");
	static PI_I = new this.CCDSC("i").scale(Math.PI);
	static #fixList(x) { //make this list contain CCDSC. mutator
		for (let i = 0; i < x.length; i++) {
			if (typeof x[i] == "number") {
				if (Number.isNaN(x[i]) || Math.abs(x[i]) == Infinity)
					throw new Error("List contains an infinite/NaN number");
				x[i] = new this.CCDSC(x[i]);
			}
			if (typeof x[i] == "string")
				try {
					x[i] = new this.CCDSC(x[i]);
				} catch(e) {
					throw new Error("Couldn't construct a CCDSC out of " + JSON.stringify(x[i]));
				}
		}
		return x;
	}
	static add(...x) {
		this.#fixList(x);
		return x.reduce((a,b)=>a.plus(b),this.ZERO);
	}
	static multiply(...x) {
		this.#fixList(x);
		return x.reduce((a,b)=>a.times(b),this.ONE);
	}
	static subtract(a, b) {
		[a, b] = this.#fixList([a, b]);
		return a.minus(b);
	}
	static divide(a, b) {
		[a, b] = this.#fixList([a, b]);
		try {
			return a.over(b);
		} catch(e) {
			throw new Error("Division by " + b.toString() + " is undefined");
		}
	}
	static abs(x) {
		[x] = this.#fixList([x]);
		return x.magnitude();
	}
	static exp(x) {
		[x] = this.#fixList([x]);
		return x.exp();
	}
	static log(x) {
		[x] = this.#fixList([x]);
		try {
			return x.ln();
		} catch(e) {
			throw new Error("Logarithm of " + x.toString() + " is undefined");
		}
	}
	static pow(a, b) {
		try {
			//a ^ b = e ^ (b ln a)
			[a] = this.#fixList([a]);
			if (typeof b == "number") {
				return a.ln().scale(b).exp();
			}
			[b] = this.#fixList([b]);
			return a.ln().times(b).exp();
			
		} catch(e) {
			//TODO: check if b is integer, do pow the old fashioned way (for cases like e^n)
			throw new Error("Couldn't calculate power because log(" + a.toString() + ") is undefined");
		}
	}
	static sign(x) {
		[x] = this.#fixList([x]);
		return x.scale(1 / x.magnitude());
	}
	static sqrt(x) {
		return this.pow(x, 0.5);
	}
	static cbrt(x) {
		return this.pow(x, 1/3);
	}
	static expm1(x) {
		[x] = this.#fixList([x]);
		return x.exp().minus(this.ONE);
	}
	static log10(x) {
		[x] = this.#fixList([x]);
		return x.ln().scale(1 / Math.LN10);
	}
	static log2(x) {
		[x] = this.#fixList([x]);
		return x.ln().scale(1 / Math.LN2);
	}
	static log1p(x) {
		[x] = this.#fixList([x]);
		return x.plus(this.ONE).ln();
	}
	static sin(x) {
		[x] = this.#fixList([x]);
		return x.times(this.I).exp().minus(x.times(this.NEGATIVE_I).exp()).times(this.NEGATIVE_HALF_I);
	}
	static cos(x) {
		[x] = this.#fixList([x]);
		return x.times(this.I).exp().plus(x.times(this.NEGATIVE_I).exp()).scale(0.5);
	}
	static sinh(x) {
		[x] = this.#fixList([x]);
		return x.exp().minus(x.scale(-1).exp()).scale(0.5);
	}
	static cosh(x) {
		[x] = this.#fixList([x]);
		return x.exp().plus(x.scale(-1).exp()).scale(0.5);
	}
	static tan(x) {
		return this.sin(x).over(this.cos(x));
	}
	static tanh(x) {
		return this.sinh(x).over(this.cosh(x));
	}
	static sec(x) {
		return this.ONE.over(this.cos(x));
	}
	static sech(x) {
		return this.ONE.over(this.cosh(x));
	}
	static csc(x) {
		return this.ONE.over(this.sin(x));
	}
	static csch(x) {
		return this.ONE.over(this.sinh(x));
	}
	static cot(x) {
		return this.cos(x).over(this.sin(x));
	}
	static coth(x) {
		return this.cosh(x).over(this.sinh(x));
	}
	static asin(x) {
		[x] = this.#fixList([x]);
		return this.sqrt(this.ONE.minus(x.times(x))).plus(this.I.times(x)).ln().times(this.NEGATIVE_I);
	}
	static asinh(x) {
		[x] = this.#fixList([x]);
		return this.sqrt(x.times(x).plus(this.ONE)).plus(x).ln();
	}
	static acos(x) {
		[x] = this.#fixList([x]);
		return this.sqrt(this.ONE.minus(x.times(x))).times(this.I).plus(x).ln().times(this.NEGATIVE_I);
	}
	static acosh(x) {
		[x] = this.#fixList([x]);
		return this.sqrt(x.times(x).minus(this.ONE)).plus(x).ln();
	}
	static atan(x) {
		[x] = this.#fixList([x]);
		// i/2 ln((i+2)/(i-2))
		return this.I.plus(x).over(this.I.minus(x)).ln().times(this.HALF_I);
	}
	static atanh(x) {
		[x] = this.#fixList([x]);
		return this.ONE.plus(x).over(this.ONE.minus(x)).ln().scale(0.5);
	}
	static asec(x) {
		[x] = this.#fixList([x]);
		return this.sqrt(this.ONE.minus(this.ONE.over(x.times(x)))).times(this.I).plus(this.ONE.over(x)).ln().times(this.NEGATIVE_I);
	}
	static asech(x) {
		[x] = this.#fixList([x]);
		return this.ONE.over(x).minus(this.sqrt(this.ONE.over(x.times(x)).minus(this.ONE))).ln();
	}
	static acsc(x) {
		[x] = this.#fixList([x]);
		return this.sqrt(this.ONE.minus(this.ONE.over(x.times(x)))).plus(this.I.over(x)).ln().times(this.NEGATIVE_I);
	}
	static acsch(x) {
		[x] = this.#fixList([x]);
		return this.ONE.over(x).plus(this.sqrt(this.ONE.over(x.times(x)).plus(this.ONE))).ln();
	}
	static acot(x) {
		[x] = this.#fixList([x]);
		return x.plus(this.I).over(x.minus(this.I)).ln().times(this.NEGATIVE_HALF_I);
	}
	static acoth(x) {
		[x] = this.#fixList([x]);
		return x.plus(this.ONE).over(x.minus(this.ONE)).ln().scale(0.5);
	}
}