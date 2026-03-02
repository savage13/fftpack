// src/fftpack.ts
function fft(x) {
  const n = x.length;
  if (n == 1) {
    return [x[0], 0];
  }
  let { wa, ifac } = cffti1(n);
  let ch = Array(2 * n).fill(0);
  let x2 = x.flatMap((v) => [v, 0]);
  return cfftf1(n, x2, ch, wa, ifac, -1);
}
function cfftf1(n, c, ch, wa, ifac, isign) {
  let ix2, ix3, ix4;
  let nf = ifac[1];
  let cinput;
  let coutput;
  let na = 0;
  let l1 = 1;
  let iw = 0;
  for (let k1 = 2; k1 <= nf + 1; k1++) {
    let ip = ifac[k1];
    let l2 = ip * l1;
    let ido = Math.floor(n / l2);
    let idot = ido + ido;
    let idl1 = idot * l1;
    if (na) {
      cinput = ch;
      coutput = c;
    } else {
      cinput = c;
      coutput = ch;
    }
    switch (ip) {
      case 4:
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        passf4(
          idot,
          l1,
          cinput,
          coutput,
          wa.slice(iw, iw + idot),
          wa.slice(ix2, ix2 + idot),
          wa.slice(ix3, ix3 + idot),
          isign
        );
        na = na == 0 ? 1 : 0;
        break;
      case 2:
        passf2(
          idot,
          l1,
          cinput,
          coutput,
          wa.slice(iw, iw + idot),
          isign
        );
        na = na == 0 ? 1 : 0;
        break;
      case 3:
        ix2 = iw + idot;
        passf3(
          idot,
          l1,
          cinput,
          coutput,
          wa.slice(iw, iw + idot),
          wa.slice(ix2, ix2 + idot),
          isign
        );
        na = na == 0 ? 1 : 0;
        break;
      case 5:
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        ix4 = ix3 + idot;
        passf5(
          idot,
          l1,
          cinput,
          coutput,
          wa.slice(iw, iw + idot),
          wa.slice(ix2, ix2 + idot),
          wa.slice(ix3, ix3 + idot),
          wa.slice(ix4, ix4 + idot),
          isign
        );
        na = na == 0 ? 1 : 0;
        break;
      default:
        let nac = passf(idot, ip, l1, idl1, cinput, coutput, wa.slice(iw), isign);
        if (nac != 0) {
          na = na == 0 ? 1 : 0;
        }
        break;
    }
    l1 = l2;
    iw += (ip - 1) * idot;
  }
  if (na == 0) {
    return c;
  }
  for (let i = 0; i < 2 * n; i++) {
    c[i] = ch[i];
  }
  return c;
}
function passf2(ido, l1, cc, ch, wa1, isign) {
  if (ido <= 2) {
    for (let k = 0; k < l1; k++) {
      let ah = k * ido;
      let ac = 2 * k * ido;
      ch[ah] = cc[ac] + cc[ac + ido];
      ch[ah + ido * l1] = cc[ac] - cc[ac + ido];
      ch[ah + 1] = cc[ac + 1] + cc[ac + ido + 1];
      ch[ah + ido * l1 + 1] = cc[ac + 1] - cc[ac + ido + 1];
    }
    return;
  } else {
    for (let k = 0; k < l1; k++) {
      for (let i = 0; i < ido - 1; i += 2) {
        let ah = i + k * ido;
        let ac = i + 2 * k * ido;
        ch[ah] = cc[ac] + cc[ac + ido];
        let tr2 = cc[ac] - cc[ac + ido];
        ch[ah + 1] = cc[ac + 1] + cc[ac + 1 + ido];
        let ti2 = cc[ac + 1] - cc[ac + 1 + ido];
        ch[ah + l1 * ido + 1] = wa1[i] * ti2 + isign * wa1[i + 1] * tr2;
        ch[ah + l1 * ido] = wa1[i] * tr2 - isign * wa1[i + 1] * ti2;
      }
    }
  }
}
function passf3(ido, l1, cc, ch, wa1, wa2, isign) {
  const taur = -0.5;
  const taui = Math.sqrt(3) / 2;
  if (ido == 2) {
    for (let k = 1; k <= l1; k++) {
      let ac = (3 * k - 2) * ido;
      let tr2 = cc[ac] + cc[ac + ido];
      let cr2 = cc[ac - ido] + taur * tr2;
      let ah = (k - 1) * ido;
      ch[ah] = cc[ac - ido] + tr2;
      let ti2 = cc[ac + 1] + cc[ac + ido + 1];
      let ci2 = cc[ac - ido + 1] + taur * ti2;
      ch[ah + 1] = cc[ac - ido + 1] + ti2;
      let cr3 = isign * taui * (cc[ac] - cc[ac + ido]);
      let ci3 = isign * taui * (cc[ac + 1] - cc[ac + ido + 1]);
      ch[ah + l1 * ido] = cr2 - ci3;
      ch[ah + 2 * l1 * ido] = cr2 + ci3;
      ch[ah + l1 * ido + 1] = ci2 + cr3;
      ch[ah + 2 * l1 * ido + 1] = ci2 - cr3;
    }
  } else {
    for (let k = 1; k <= l1; k++) {
      for (let i = 0; i < ido - 1; i += 2) {
        let ac = i + (3 * k - 2) * ido;
        let tr2 = cc[ac] + cc[ac + ido];
        let cr2 = cc[ac - ido] + taur * tr2;
        let ah = i + (k - 1) * ido;
        ch[ah] = cc[ac - ido] + tr2;
        let ti2 = cc[ac + 1] + cc[ac + ido + 1];
        let ci2 = cc[ac - ido + 1] + taur * ti2;
        ch[ah + 1] = cc[ac - ido + 1] + ti2;
        let cr3 = isign * taui * (cc[ac] - cc[ac + ido]);
        let ci3 = isign * taui * (cc[ac + 1] - cc[ac + ido + 1]);
        let dr2 = cr2 - ci3;
        let dr3 = cr2 + ci3;
        let di2 = ci2 + cr3;
        let di3 = ci2 - cr3;
        ch[ah + l1 * ido + 1] = wa1[i] * di2 + isign * wa1[i + 1] * dr2;
        ch[ah + l1 * ido] = wa1[i] * dr2 - isign * wa1[i + 1] * di2;
        ch[ah + 2 * l1 * ido + 1] = wa2[i] * di3 + isign * wa2[i + 1] * dr3;
        ch[ah + 2 * l1 * ido] = wa2[i] * dr3 - isign * wa2[i + 1] * di3;
      }
    }
  }
}
function passf4(ido, l1, cc, ch, wa1, wa2, wa3, isign) {
  if (ido == 2) {
    for (let k = 0; k < l1; k++) {
      let ac = 4 * k * ido + 1;
      let ti1 = cc[ac] - cc[ac + 2 * ido];
      let ti2 = cc[ac] + cc[ac + 2 * ido];
      let tr4 = cc[ac + 3 * ido] - cc[ac + ido];
      let ti3 = cc[ac + ido] + cc[ac + 3 * ido];
      let tr1 = cc[ac - 1] - cc[ac + 2 * ido - 1];
      let tr2 = cc[ac - 1] + cc[ac + 2 * ido - 1];
      let ti4 = cc[ac + ido - 1] - cc[ac + 3 * ido - 1];
      let tr3 = cc[ac + ido - 1] + cc[ac + 3 * ido - 1];
      let ah = k * ido;
      ch[ah] = tr2 + tr3;
      ch[ah + 2 * l1 * ido] = tr2 - tr3;
      ch[ah + 1] = ti2 + ti3;
      ch[ah + 2 * l1 * ido + 1] = ti2 - ti3;
      ch[ah + l1 * ido] = tr1 + isign * tr4;
      ch[ah + 3 * l1 * ido] = tr1 - isign * tr4;
      ch[ah + l1 * ido + 1] = ti1 + isign * ti4;
      ch[ah + 3 * l1 * ido + 1] = ti1 - isign * ti4;
    }
  } else {
    for (let k = 0; k < l1; k++) {
      for (let i = 0; i < ido - 1; i += 2) {
        let ac = i + 1 + 4 * k * ido;
        let ti1 = cc[ac] - cc[ac + 2 * ido];
        let ti2 = cc[ac] + cc[ac + 2 * ido];
        let ti3 = cc[ac + ido] + cc[ac + 3 * ido];
        let tr4 = cc[ac + 3 * ido] - cc[ac + ido];
        let tr1 = cc[ac - 1] - cc[ac + 2 * ido - 1];
        let tr2 = cc[ac - 1] + cc[ac + 2 * ido - 1];
        let ti4 = cc[ac + ido - 1] - cc[ac + 3 * ido - 1];
        let tr3 = cc[ac + ido - 1] + cc[ac + 3 * ido - 1];
        let ah = i + k * ido;
        ch[ah] = tr2 + tr3;
        let cr3 = tr2 - tr3;
        ch[ah + 1] = ti2 + ti3;
        let ci3 = ti2 - ti3;
        let cr2 = tr1 + isign * tr4;
        let cr4 = tr1 - isign * tr4;
        let ci2 = ti1 + isign * ti4;
        let ci4 = ti1 - isign * ti4;
        ch[ah + l1 * ido] = wa1[i] * cr2 - isign * wa1[i + 1] * ci2;
        ch[ah + l1 * ido + 1] = wa1[i] * ci2 + isign * wa1[i + 1] * cr2;
        ch[ah + 2 * l1 * ido] = wa2[i] * cr3 - isign * wa2[i + 1] * ci3;
        ch[ah + 2 * l1 * ido + 1] = wa2[i] * ci3 + isign * wa2[i + 1] * cr3;
        ch[ah + 3 * l1 * ido] = wa3[i] * cr4 - isign * wa3[i + 1] * ci4;
        ch[ah + 3 * l1 * ido + 1] = wa3[i] * ci4 + isign * wa3[i + 1] * cr4;
      }
    }
  }
}
function passf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4, isign) {
  const tr11 = 1 / (Math.sqrt(5) + 1);
  const ti11 = 0.5 * Math.sqrt(0.5 * (5 + Math.sqrt(5)));
  const tr12 = -1 / (Math.sqrt(5) - 1);
  const ti12 = Math.sqrt(5 / (2 * (5 + Math.sqrt(5))));
  if (ido == 2) {
    for (let k = 1; k <= l1; ++k) {
      let ac = (5 * k - 4) * ido + 1;
      let ti5 = cc[ac] - cc[ac + 3 * ido];
      let ti2 = cc[ac] + cc[ac + 3 * ido];
      let ti4 = cc[ac + ido] - cc[ac + 2 * ido];
      let ti3 = cc[ac + ido] + cc[ac + 2 * ido];
      let tr5 = cc[ac - 1] - cc[ac + 3 * ido - 1];
      let tr2 = cc[ac - 1] + cc[ac + 3 * ido - 1];
      let tr4 = cc[ac + ido - 1] - cc[ac + 2 * ido - 1];
      let tr3 = cc[ac + ido - 1] + cc[ac + 2 * ido - 1];
      let ah = (k - 1) * ido;
      ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
      ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
      let cr2 = cc[ac - ido - 1] + tr11 * tr2 + tr12 * tr3;
      let ci2 = cc[ac - ido] + tr11 * ti2 + tr12 * ti3;
      let cr3 = cc[ac - ido - 1] + tr12 * tr2 + tr11 * tr3;
      let ci3 = cc[ac - ido] + tr12 * ti2 + tr11 * ti3;
      let cr5 = isign * (ti11 * tr5 + ti12 * tr4);
      let ci5 = isign * (ti11 * ti5 + ti12 * ti4);
      let cr4 = isign * (ti12 * tr5 - ti11 * tr4);
      let ci4 = isign * (ti12 * ti5 - ti11 * ti4);
      ch[ah + l1 * ido] = cr2 - ci5;
      ch[ah + 4 * l1 * ido] = cr2 + ci5;
      ch[ah + l1 * ido + 1] = ci2 + cr5;
      ch[ah + 2 * l1 * ido + 1] = ci3 + cr4;
      ch[ah + 2 * l1 * ido] = cr3 - ci4;
      ch[ah + 3 * l1 * ido] = cr3 + ci4;
      ch[ah + 3 * l1 * ido + 1] = ci3 - cr4;
      ch[ah + 4 * l1 * ido + 1] = ci2 - cr5;
    }
  } else {
    for (let k = 1; k <= l1; k++) {
      for (let i = 0; i < ido - 1; i += 2) {
        let ac = i + 1 + (k * 5 - 4) * ido;
        let ti5 = cc[ac] - cc[ac + 3 * ido];
        let ti2 = cc[ac] + cc[ac + 3 * ido];
        let ti4 = cc[ac + ido] - cc[ac + 2 * ido];
        let ti3 = cc[ac + ido] + cc[ac + 2 * ido];
        let tr5 = cc[ac - 1] - cc[ac + 3 * ido - 1];
        let tr2 = cc[ac - 1] + cc[ac + 3 * ido - 1];
        let tr4 = cc[ac + ido - 1] - cc[ac + 2 * ido - 1];
        let tr3 = cc[ac + ido - 1] + cc[ac + 2 * ido - 1];
        let ah = i + (k - 1) * ido;
        ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
        ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
        let cr2 = cc[ac - ido - 1] + tr11 * tr2 + tr12 * tr3;
        let ci2 = cc[ac - ido] + tr11 * ti2 + tr12 * ti3;
        let cr3 = cc[ac - ido - 1] + tr12 * tr2 + tr11 * tr3;
        let ci3 = cc[ac - ido] + tr12 * ti2 + tr11 * ti3;
        let cr5 = isign * (ti11 * tr5 + ti12 * tr4);
        let ci5 = isign * (ti11 * ti5 + ti12 * ti4);
        let cr4 = isign * (ti12 * tr5 - ti11 * tr4);
        let ci4 = isign * (ti12 * ti5 - ti11 * ti4);
        let dr3 = cr3 - ci4;
        let dr4 = cr3 + ci4;
        let di3 = ci3 + cr4;
        let di4 = ci3 - cr4;
        let dr5 = cr2 + ci5;
        let dr2 = cr2 - ci5;
        let di5 = ci2 - cr5;
        let di2 = ci2 + cr5;
        ch[ah + l1 * ido] = wa1[i] * dr2 - isign * wa1[i + 1] * di2;
        ch[ah + l1 * ido + 1] = wa1[i] * di2 + isign * wa1[i + 1] * dr2;
        ch[ah + 2 * l1 * ido] = wa2[i] * dr3 - isign * wa2[i + 1] * di3;
        ch[ah + 2 * l1 * ido + 1] = wa2[i] * di3 + isign * wa2[i + 1] * dr3;
        ch[ah + 3 * l1 * ido] = wa3[i] * dr4 - isign * wa3[i + 1] * di4;
        ch[ah + 3 * l1 * ido + 1] = wa3[i] * di4 + isign * wa3[i + 1] * dr4;
        ch[ah + 4 * l1 * ido] = wa4[i] * dr5 - isign * wa4[i + 1] * di5;
        ch[ah + 4 * l1 * ido + 1] = wa4[i] * di5 + isign * wa4[i + 1] * dr5;
      }
    }
  }
}
function passf(ido, ip, l1, idl1, cc, ch, wa, isign) {
  let nac;
  let idot = ido / 2;
  let ipph = (ip + 1) / 2;
  let idp = ip * ido;
  if (ido >= l1) {
    for (let j = 1; j < ipph; j++) {
      let jc = ip - j;
      for (let k = 0; k < l1; k++) {
        for (let i = 0; i < ido; i++) {
          ch[i + (k + j * l1) * ido] = cc[i + (j + k * ip) * ido] + cc[i + (jc + k * ip) * ido];
          ch[i + (k + jc * l1) * ido] = cc[i + (j + k * ip) * ido] - cc[i + (jc + k * ip) * ido];
        }
      }
    }
    for (let k = 0; k < l1; k++) {
      for (let i = 0; i < ido; i++) {
        ch[i + k * ido] = cc[i + k * ip * ido];
      }
    }
  } else {
    for (let j = 1; j < ipph; j++) {
      let jc = ip - j;
      for (let i = 0; i < ido; i++) {
        for (let k = 0; k < l1; k++) {
          ch[i + (k + j * l1) * ido] = cc[i + (j + k * ip) * ido] + cc[i + (jc + k * ip) * ido];
          ch[i + (k + jc * l1) * ido] = cc[i + (j + k * ip) * ido] - cc[i + (jc + k * ip) * ido];
        }
      }
    }
    for (let i = 0; i < ido; i++) {
      for (let k = 0; k < l1; k++) {
        ch[i + k * ido] = cc[i + k * ip * ido];
      }
    }
  }
  let idl = 2 - ido;
  let inc = 0;
  for (let l = 1; l < ipph; l++) {
    let lc = ip - l;
    idl += ido;
    for (let ik = 0; ik < idl1; ik++) {
      cc[ik + l * idl1] = ch[ik] + wa[idl - 2] * ch[ik + idl1];
      cc[ik + lc * idl1] = isign * wa[idl - 1] * ch[ik + (ip - 1) * idl1];
    }
    let idlj = idl;
    inc += ido;
    for (let j = 2; j < ipph; j++) {
      let jc = ip - j;
      idlj += inc;
      if (idlj > idp) idlj -= idp;
      let war = wa[idlj - 2];
      let wai = wa[idlj - 1];
      for (let ik = 0; ik < idl1; ik++) {
        cc[ik + l * idl1] += war * ch[ik + j * idl1];
        cc[ik + lc * idl1] += isign * wai * ch[ik + jc * idl1];
      }
    }
  }
  for (let j = 1; j < ipph; j++) {
    for (let ik = 0; ik < idl1; ik++) {
      ch[ik] += ch[ik + j * idl1];
    }
  }
  for (let j = 1; j < ipph; j++) {
    let jc = ip - j;
    for (let ik = 1; ik < idl1; ik += 2) {
      ch[ik - 1 + j * idl1] = cc[ik - 1 + j * idl1] - cc[ik + jc * idl1];
      ch[ik - 1 + jc * idl1] = cc[ik - 1 + j * idl1] + cc[ik + jc * idl1];
      ch[ik + j * idl1] = cc[ik + j * idl1] + cc[ik - 1 + jc * idl1];
      ch[ik + jc * idl1] = cc[ik + j * idl1] - cc[ik - 1 + jc * idl1];
    }
  }
  nac = 1;
  if (ido == 2) {
    return nac;
  }
  nac = 0;
  for (let ik = 0; ik < idl1; ik++) {
    cc[ik] = ch[ik];
  }
  for (let j = 1; j < ip; j++) {
    for (let k = 0; k < l1; k++) {
      cc[(k + j * l1) * ido + 0] = ch[(k + j * l1) * ido + 0];
      cc[(k + j * l1) * ido + 1] = ch[(k + j * l1) * ido + 1];
    }
  }
  if (idot <= l1) {
    let idij = 0;
    for (let j = 1; j < ip; j++) {
      idij += 2;
      for (let i = 3; i < ido; i += 2) {
        idij += 2;
        for (let k = 0; k < l1; k++) {
          cc[i - 1 + (k + j * l1) * ido] = wa[idij - 2] * ch[i - 1 + (k + j * l1) * ido] - isign * wa[idij - 1] * ch[i + (k + j * l1) * ido];
          cc[i + (k + j * l1) * ido] = wa[idij - 2] * ch[i + (k + j * l1) * ido] + isign * wa[idij - 1] * ch[i - 1 + (k + j * l1) * ido];
        }
      }
    }
  } else {
    let idj = 2 - ido;
    for (let j = 1; j < ip; j++) {
      idj += ido;
      for (let k = 0; k < l1; k++) {
        let idij = idj;
        for (let i = 3; i < ido; i += 2) {
          idij += 2;
          cc[i - 1 + (k + j * l1) * ido] = wa[idij - 2] * ch[i - 1 + (k + j * l1) * ido] - isign * wa[idij - 1] * ch[i + (k + j * l1) * ido];
          cc[i + (k + j * l1) * ido] = wa[idij - 2] * ch[i + (k + j * l1) * ido] + isign * wa[idij - 1] * ch[i - 1 + (k + j * l1) * ido];
        }
      }
    }
  }
  return nac;
}
function cffti(n) {
  if (n == 1)
    return;
  return cffti1(n);
}
function cffti1(n) {
  let wa = [];
  let ifac = factorize(n);
  let nf = ifac[1];
  const argh = 2 * Math.PI / n;
  let i = 1;
  let l1 = 1;
  for (let k1 = 1; k1 <= nf; k1++) {
    let ip = ifac[k1 + 1];
    let ld = 0;
    let l2 = l1 * ip;
    let ido = Math.floor(n / l2);
    let idot = ido + ido + 2;
    let ipm = ip - 1;
    for (let j = 1; j <= ipm; j++) {
      let i1 = i;
      wa[i - 1] = 1;
      wa[i] = 0;
      ld += l1;
      let fi = 0;
      let argld = ld * argh;
      for (let ii = 4; ii <= idot; ii += 2) {
        i += 2;
        fi += 1;
        let arg = fi * argld;
        wa[i - 1] = Math.cos(arg);
        wa[i] = Math.sin(arg);
      }
      if (ip > 5) {
        wa[i1 - 1] = wa[i - 1];
        wa[i1] = wa[i];
      }
    }
    l1 = l2;
  }
  return { wa, ifac };
}
function factorize(n) {
  let ifac = [];
  let ntryh = [3, 4, 2, 5];
  let ntry = 3;
  let nl = n;
  let nf = 0;
  let j = 0;
  let startloop = true;
  while (startloop) {
    startloop = false;
    if (j < ntryh.length) {
      ntry = ntryh[j];
    } else {
      ntry += 2;
    }
    j += 1;
    do {
      let nq = Math.floor(nl / ntry);
      let nr = nl - ntry * nq;
      if (nr != 0) {
        startloop = true;
        break;
      }
      nf += 1;
      ifac[nf + 1] = ntry;
      nl = nq;
      if (ntry == 2 && nf != 1) {
        for (let i = 2; i <= nf; i++) {
          let ib = nf - i + 2;
          ifac[ib + 1] = ifac[ib];
        }
        ifac[2] = 2;
      }
    } while (nl != 1);
  }
  ifac[0] = n;
  ifac[1] = nf;
  return ifac;
}
export {
  fft
};
//# sourceMappingURL=fftpack.js.map
