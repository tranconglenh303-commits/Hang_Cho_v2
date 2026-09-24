import { KendallModelType } from '../types/queue';

export function factorial(n: number): number {
  if (n <= 1) return 1;
  let res = 1;
  for (let i = 2; i <= n; i++) res *= i;
  return res;
}

export function calculateMM1(lambda: number, mu: number) {
  const rho = lambda / mu;
  const isStable = rho < 1;

  if (!isStable) {
    return {
      rho,
      isStable: false,
      P0: 0,
      L: Infinity,
      Lq: Infinity,
      W: Infinity,
      Wq: Infinity,
      Pw: 1,
    };
  }

  const P0 = 1 - rho;
  const L = rho / (1 - rho);
  const Lq = (rho * rho) / (1 - rho);
  const W = 1 / (mu - lambda);
  const Wq = rho / (mu - lambda);
  const Pw = rho; // PASTA: P(wait) = P(n >= 1) = rho

  return { rho, isStable: true, P0, L, Lq, W, Wq, Pw };
}

export function calculateMMC(lambda: number, mu: number, c: number) {
  const r = lambda / mu;
  const rho = r / c;
  const isStable = rho < 1;

  if (!isStable) {
    return {
      rho,
      isStable: false,
      P0: 0,
      L: Infinity,
      Lq: Infinity,
      W: Infinity,
      Wq: Infinity,
      Pw: 1,
    };
  }

  // Calculate P0
  let sum = 0;
  for (let n = 0; n < c; n++) {
    sum += Math.pow(r, n) / factorial(n);
  }
  sum += (Math.pow(r, c) / (factorial(c) * (1 - rho)));
  const P0 = 1 / sum;

  // Erlang C (Probability of waiting in queue: Pw)
  const Pw = (Math.pow(r, c) / (factorial(c) * (1 - rho))) * P0;

  // Average queue length Lq
  const Lq = (Pw * rho) / (1 - rho);
  const Wq = Lq / lambda;
  const W = Wq + 1 / mu;
  const L = lambda * W; // L = Lq + r

  return { rho, isStable: true, P0, L, Lq, W, Wq, Pw };
}

export function calculateMM1K(lambda: number, mu: number, K: number) {
  const rho = lambda / mu;
  let P0 = 0;
  let PK = 0;

  if (Math.abs(rho - 1) < 1e-6) {
    P0 = 1 / (K + 1);
    PK = 1 / (K + 1);
  } else {
    P0 = (1 - rho) / (1 - Math.pow(rho, K + 1));
    PK = (Math.pow(rho, K) * (1 - rho)) / (1 - Math.pow(rho, K + 1));
  }

  const lambdaEff = lambda * (1 - PK);
  let L = 0;
  if (Math.abs(rho - 1) < 1e-6) {
    L = K / 2;
  } else {
    L = (rho / (1 - rho)) - (((K + 1) * Math.pow(rho, K + 1)) / (1 - Math.pow(rho, K + 1)));
  }

  const W = L / lambdaEff;
  const Wq = W - 1 / mu;
  const Lq = lambdaEff * Wq;
  const Pw = 1 - P0;

  return { rho, isStable: true, P0, PK, lambdaEff, L, Lq, W, Wq, Pw };
}

export function calculateMD1(lambda: number, mu: number) {
  const rho = lambda / mu;
  const isStable = rho < 1;

  if (!isStable) {
    return {
      rho,
      isStable: false,
      P0: 0,
      L: Infinity,
      Lq: Infinity,
      W: Infinity,
      Wq: Infinity,
      Pw: 1,
    };
  }

  const P0 = 1 - rho;
  const Lq = (rho * rho) / (2 * (1 - rho));
  const Wq = Lq / lambda;
  const W = Wq + 1 / mu;
  const L = lambda * W;
  const Pw = rho;

  return { rho, isStable: true, P0, L, Lq, W, Wq, Pw };
}

export function calculateMGkLoss(lambda: number, mu: number, k: number) {
  const a = lambda / mu; // offered load
  const rho = a / k; // traffic intensity per server

  // Compute Erlang-B recursion: B(n, a) = (a * B(n-1, a)) / (n + a * B(n-1, a))
  let B = 1;
  for (let n = 1; n <= k; n++) {
    B = (a * B) / (n + a * B);
  }
  const Pk = B; // Loss probability (Erlang B)
  const Pw = Pk; // Probability of arrival being blocked

  // Compute sum_{n=0}^k a^n / n! to get P0
  let sumTerm = 1;
  let term = 1;
  for (let n = 1; n <= k; n++) {
    term = term * (a / n);
    sumTerm += term;
  }
  const P0 = 1 / sumTerm;

  const lambdaEff = lambda * (1 - Pk);
  const L = a * (1 - Pk); // avg customers in system = avg busy servers
  const Lq = 0; // No queue capacity
  const Wq = 0; // No waiting time
  const W = 1 / mu; // Service time for admitted customers

  return {
    a,
    rho,
    isStable: true, // Loss systems are always stable for any finite lambda
    Pk,
    P0,
    L,
    Lq,
    W,
    Wq,
    Pw,
    lambdaEff,
  };
}

export function calculateMG1(lambda: number, mu: number, sigma: number) {
  const rho = lambda / mu;
  const isStable = rho < 1;

  if (!isStable) {
    return {
      rho,
      isStable: false,
      P0: 0,
      L: Infinity,
      Lq: Infinity,
      W: Infinity,
      Wq: Infinity,
      Pw: 1,
    };
  }

  const P0 = 1 - rho;
  const varS = sigma * sigma;
  const Lq = (Math.pow(lambda, 2) * varS + Math.pow(rho, 2)) / (2 * (1 - rho));
  const Wq = Lq / lambda;
  const W = Wq + 1 / mu;
  const L = lambda * W;
  const Pw = rho;

  return { rho, isStable: true, P0, L, Lq, W, Wq, Pw };
}

export interface CalculationStep {
  stepNumber: number;
  stageName: string;
  description: string;
  latexFormula: string;
  numericalResult: string;
  explanation: string;
}

export function getStepByStepCalculation(
  model: KendallModelType,
  lambda: number,
  mu: number,
  k: number = 2,
  sigma: number = 0.05
): CalculationStep[] {
  if (model === 'M/M/1') {
    const res = calculateMM1(lambda, mu);
    if (!res.isStable) {
      return [
        {
          stepNumber: 1,
          stageName: 'Phép tính 1 (P₀): Xác suất quầy rảnh rỗi',
          description: 'Hệ số sử dụng ρ = λ/μ ≥ 1. Hệ thống không ổn định.',
          latexFormula: 'P_0 = 1 - \\rho = 0',
          numericalResult: `\\rho = \\frac{${lambda}}{${mu}} = ${res.rho.toFixed(3)} \\ge 1 \\implies P_0 = 0`,
          explanation: 'Hệ thống bùng nổ vô hạn vì tốc độ khách đến lớn hơn hoặc bằng tốc độ phục vụ.',
        },
        {
          stepNumber: 2,
          stageName: 'Phép tính 2 (Lq): Số lượng khách trung bình trong hàng chờ',
          description: 'Hàng chờ kéo dài vô hạn khi ρ ≥ 1.',
          latexFormula: 'L_q = \\frac{\\rho^2}{1 - \\rho} \\to \\infty',
          numericalResult: 'L_q = \\infty \\text{ (Vô hạn)}',
          explanation: 'Hàng chờ không bao giờ giải tỏa hết.',
        },
        {
          stepNumber: 3,
          stageName: 'Phép tính 3 (L): Số lượng khách trung bình trong toàn hệ thống',
          description: 'Số khách trong hệ thống tiến tới vô cùng.',
          latexFormula: 'L = \\frac{\\rho}{1 - \\rho} \\to \\infty',
          numericalResult: 'L = \\infty \\text{ (Vô hạn)}',
          explanation: 'Khách đến dồn ứ vô hạn.',
        },
        {
          stepNumber: 4,
          stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình trong hàng chờ',
          description: 'Thời gian chờ đợi tăng vô hạn.',
          latexFormula: 'W_q = \\frac{L_q}{\\lambda} \\to \\infty',
          numericalResult: 'W_q = \\infty \\text{ (Vô hạn)}',
          explanation: 'Khách hàng phải chờ vô tận.',
        },
        {
          stepNumber: 5,
          stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình trong hệ thống',
          description: 'Thời gian lưu lại tăng vô hạn.',
          latexFormula: 'W = \\frac{L}{\\lambda} \\to \\infty',
          numericalResult: 'W = \\infty \\text{ (Vô hạn)}',
          explanation: 'Hệ thống tắc nghẽn hoàn toàn.',
        },
        {
          stepNumber: 6,
          stageName: 'Phép tính 6 (Pw): Xác suất khách đến phải chờ',
          description: 'Hệ thống luôn đầy.',
          latexFormula: 'P_w = 1 - P_0 = 100\\%',
          numericalResult: 'P_w = 100\\%',
          explanation: '100% khách đến đều phải chờ.',
        },
      ];
    }

    return [
      {
        stepNumber: 1,
        stageName: 'Phép tính 1 (P₀): Xác suất không có khách trong hệ thống (Hệ thống rảnh)',
        description: 'Xác suất để hệ thống không có bất kỳ khách hàng nào, quầy phục vụ hoàn toàn rảnh rỗi.',
        latexFormula: 'P_0 = 1 - \\rho = 1 - \\frac{\\lambda}{\\mu}',
        numericalResult: `P_0 = 1 - \\frac{${lambda}}{${mu}} = 1 - ${res.rho.toFixed(3)} = ${res.P0.toFixed(4)} \\quad (${(res.P0 * 100).toFixed(2)}\\%)`,
        explanation: `Quầy rảnh rỗi ${(res.P0 * 100).toFixed(2)}% thời gian hoạt động. Quầy đang bận ${(res.rho * 100).toFixed(2)}% thời gian.`,
      },
      {
        stepNumber: 2,
        stageName: 'Phép tính 2 (Lq): Số lượng khách trung bình trong hàng chờ',
        description: 'Số lượng khách hàng trung bình đang đứng đợi trong hàng (chưa được phục vụ).',
        latexFormula: 'L_q = \\frac{\\rho^2}{1 - \\rho} = \\frac{\\lambda^2}{\\mu(\\mu - \\lambda)}',
        numericalResult: `L_q = \\frac{${res.rho.toFixed(3)}^2}{1 - ${res.rho.toFixed(3)}} = ${res.Lq.toFixed(4)} \\text{ khách}`,
        explanation: `Trung bình tại bất kỳ thời điểm nào có ${res.Lq.toFixed(2)} khách đang xếp hàng đợi tới lượt.`,
      },
      {
        stepNumber: 3,
        stageName: 'Phép tính 3 (L): Số lượng khách trung bình trong toàn hệ thống',
        description: 'Tổng số khách đang có mặt trong hệ thống (bao gồm cả khách đang đợi và khách đang được phục vụ tại quầy).',
        latexFormula: 'L = L_q + \\rho = \\frac{\\rho}{1 - \\rho} = \\frac{\\lambda}{\\mu - \\lambda}',
        numericalResult: `L = \\frac{${res.rho.toFixed(3)}}{1 - ${res.rho.toFixed(3)}} = ${res.L.toFixed(4)} \\text{ khách}`,
        explanation: `Trung bình có ${res.L.toFixed(2)} khách trong toàn bộ cửa hàng (${res.Lq.toFixed(2)} đang đợi + ${res.rho.toFixed(2)} đang giao dịch).`,
      },
      {
        stepNumber: 4,
        stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình trong hàng chờ',
        description: 'Thời gian một khách hàng phải xếp hàng đợi trước khi bắt đầu được nhân viên phục vụ (Định luật Little: Wq = Lq / λ).',
        latexFormula: 'W_q = \\frac{L_q}{\\lambda} = \\frac{\\rho}{\\mu - \\lambda} = \\frac{\\lambda}{\\mu(\\mu - \\lambda)}',
        numericalResult: `W_q = \\frac{${res.Lq.toFixed(4)}}{${lambda}} = ${res.Wq.toFixed(4)} \\text{ giờ} = ${(res.Wq * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Mỗi khách hàng đến phải chờ trung bình ${(res.Wq * 60).toFixed(2)} phút trước khi được nhân viên gọi tên.`,
      },
      {
        stepNumber: 5,
        stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình trong toàn hệ thống',
        description: 'Tổng thời gian một khách hàng trải qua từ khi bước vào hệ thống cho tới khi hoàn tất phục vụ ra về (W = Wq + 1/μ).',
        latexFormula: 'W = W_q + \\frac{1}{\\mu} = \\frac{1}{\\mu - \\lambda} = \\frac{L}{\\lambda}',
        numericalResult: `W = ${res.Wq.toFixed(4)} + \\frac{1}{${mu}} = ${res.W.toFixed(4)} \\text{ giờ} = ${(res.W * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Tổng thời gian từ lúc đến đến khi rời đi là ${(res.W * 60).toFixed(2)} phút (bao gồm ${(res.Wq * 60).toFixed(2)} phút chờ và ${((1 / mu) * 60).toFixed(2)} phút phục vụ).`,
      },
      {
        stepNumber: 6,
        stageName: 'Phép tính 6 (Pw): Xác suất khách hàng khi đến phải chờ',
        description: 'Xác suất một khách hàng mới đến thấy quầy đang bận phục vụ người khác và buộc phải vào hàng chờ.',
        latexFormula: 'P_w = P(n \\ge 1) = 1 - P_0 = \\rho = \\frac{\\lambda}{\\mu}',
        numericalResult: `P_w = 1 - ${res.P0.toFixed(4)} = ${res.Pw.toFixed(4)} \\quad (${(res.Pw * 100).toFixed(2)}\\%)`,
        explanation: `Có ${(res.Pw * 100).toFixed(2)}% khả năng một khách hàng khi bước vào sẽ phải xếp hàng chờ đợi.`,
      },
    ];
  }

  if (model === 'M/M/k') {
    const res = calculateMMC(lambda, mu, k);
    const r = lambda / mu;

    if (!res.isStable) {
      return [
        {
          stepNumber: 1,
          stageName: 'Phép tính 1 (P₀): Xác suất hệ thống rảnh',
          description: 'Hệ số sử dụng ρ = λ/(k·μ) ≥ 1. Hệ thống không ổn định.',
          latexFormula: 'P_0 = 0',
          numericalResult: `\\rho = \\frac{${lambda}}{${k} \\times ${mu}} = ${res.rho.toFixed(3)} \\ge 1`,
          explanation: 'Hệ thống bùng nổ vô hạn vì tốc độ khách đến lớn hơn hoặc bằng tổng công suất của k quầy.',
        },
        {
          stepNumber: 2,
          stageName: 'Phép tính 2 (Lq): Chiều dài hàng chờ',
          description: 'Hàng chờ kéo dài vô hạn.',
          latexFormula: 'L_q \\to \\infty',
          numericalResult: 'L_q = \\infty \\text{ (Vô hạn)}',
          explanation: 'Hàng chờ dồn ứ vô hạn.',
        },
        {
          stepNumber: 3,
          stageName: 'Phép tính 3 (L): Số khách trong toàn hệ thống',
          description: 'Hệ thống quá tải.',
          latexFormula: 'L \\to \\infty',
          numericalResult: 'L = \\infty',
          explanation: 'Quá tải toàn diện.',
        },
        {
          stepNumber: 4,
          stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình',
          description: 'Thời gian chờ vô tận.',
          latexFormula: 'W_q \\to \\infty',
          numericalResult: 'W_q = \\infty',
          explanation: 'Khách hàng không bao giờ được phục vụ hết.',
        },
        {
          stepNumber: 5,
          stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình',
          description: 'Thời gian lưu lại vô tận.',
          latexFormula: 'W \\to \\infty',
          numericalResult: 'W = \\infty',
          explanation: 'Hệ thống sụp đổ.',
        },
        {
          stepNumber: 6,
          stageName: 'Phép tính 6 (Pw): Xác suất phải xếp hàng',
          description: 'Tất cả quầy đều bận 100%.',
          latexFormula: 'P_w = 100\\%',
          numericalResult: 'P_w = 100\\%',
          explanation: 'Mọi khách đến đều phải chờ.',
        },
      ];
    }

    return [
      {
        stepNumber: 1,
        stageName: 'Phép tính 1 (P₀): Xác suất không có khách trong hệ thống (Tất cả quầy đều rảnh)',
        description: `Mô hình M/M/${k}: Xác suất toàn bộ ${k} quầy phục vụ đều không có khách nào (cả hệ thống rảnh).`,
        latexFormula: 'P_0 = \\left[ \\sum_{n=0}^{k-1} \\frac{r^n}{n!} + \\frac{r^k}{k!(1 - \\rho)} \\right]^{-1} \\quad \\left( r = \\frac{\\lambda}{\\mu} \\right)',
        numericalResult: `r = \\frac{${lambda}}{${mu}} = ${r.toFixed(2)}, \\; \\rho = \\frac{${r.toFixed(2)}}{${k}} = ${res.rho.toFixed(3)} \\implies P_0 = ${res.P0.toFixed(4)} \\quad (${(res.P0 * 100).toFixed(2)}\\%)`,
        explanation: `Xác suất để tất cả ${k} quầy cùng rảnh rỗi đồng thời là ${(res.P0 * 100).toFixed(2)}%.`,
      },
      {
        stepNumber: 2,
        stageName: 'Phép tính 2 (Lq): Số lượng khách trung bình trong hàng chờ',
        description: 'Số lượng khách hàng trung bình đứng đợi trong hàng chờ chung (Erlang C).',
        latexFormula: 'L_q = \\frac{P_w \\cdot \\rho}{1 - \\rho} = \\frac{P_0 \\cdot r^k \\cdot \\rho}{k!(1 - \\rho)^2}',
        numericalResult: `L_q = \\frac{${res.Pw.toFixed(4)} \\times ${res.rho.toFixed(3)}}{1 - ${res.rho.toFixed(3)}} = ${res.Lq.toFixed(4)} \\text{ khách}`,
        explanation: `Trung bình có ${res.Lq.toFixed(2)} khách đang xếp hàng đợi ở hàng chờ dùng chung.`,
      },
      {
        stepNumber: 3,
        stageName: 'Phép tính 3 (L): Số lượng khách trung bình trong toàn hệ thống',
        description: `Tổng số khách trong hệ thống gồm khách đợi trong hàng và khách đang được phục vụ tại ${k} quầy.`,
        latexFormula: 'L = L_q + r = L_q + \\frac{\\lambda}{\\mu} = \\lambda W',
        numericalResult: `L = ${res.Lq.toFixed(4)} + ${r.toFixed(4)} = ${res.L.toFixed(4)} \\text{ khách}`,
        explanation: `Trung bình có ${res.L.toFixed(2)} khách trong toàn bộ hệ thống (${res.Lq.toFixed(2)} khách đợi + ${r.toFixed(2)} khách đang giao dịch tại các quầy).`,
      },
      {
        stepNumber: 4,
        stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình trong hàng chờ',
        description: 'Thời gian một khách hàng trung bình phải đợi trong hàng trước khi có 1 trong các quầy trống gọi vào.',
        latexFormula: 'W_q = \\frac{L_q}{\\lambda} = \\frac{P_w}{k\\mu - \\lambda}',
        numericalResult: `W_q = \\frac{${res.Lq.toFixed(4)}}{${lambda}} = ${res.Wq.toFixed(4)} \\text{ giờ} = ${(res.Wq * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Khách hàng chỉ phải đợi trung bình ${(res.Wq * 60).toFixed(2)} phút trong hàng chờ chung.`,
      },
      {
        stepNumber: 5,
        stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình trong toàn hệ thống',
        description: 'Tổng thời gian một khách hàng lưu lại từ khi đến cho tới khi phục vụ xong.',
        latexFormula: 'W = W_q + \\frac{1}{\\mu} = \\frac{L}{\\lambda}',
        numericalResult: `W = ${res.Wq.toFixed(4)} + \\frac{1}{${mu}} = ${res.W.toFixed(4)} \\text{ giờ} = ${(res.W * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Tổng thời gian khách hàng lưu lại là ${(res.W * 60).toFixed(2)} phút (gồm ${(res.Wq * 60).toFixed(2)} phút chờ và ${((1 / mu) * 60).toFixed(2)} phút phục vụ).`,
      },
      {
        stepNumber: 6,
        stageName: 'Phép tính 6 (Pw): Xác suất khách hàng khi đến phải chờ (Công thức Erlang C)',
        description: `Xác suất một khách hàng khi đến thấy toàn bộ ${k} quầy đều đang bận và buộc phải xếp hàng.`,
        latexFormula: 'P_w = C(k, r) = \\frac{\\frac{r^k}{k!(1 - \\rho)}}{\\sum_{n=0}^{k-1} \\frac{r^n}{n!} + \\frac{r^k}{k!(1 - \\rho)}} = \\frac{r^k}{k!(1 - \\rho)} \\cdot P_0',
        numericalResult: `P_w = ${res.Pw.toFixed(4)} \\quad (${(res.Pw * 100).toFixed(2)}\\%)`,
        explanation: `Có ${(res.Pw * 100).toFixed(2)}% khả năng một khách hàng khi đến phải đứng đợi vì cả ${k} quầy đều đang bận.`,
      },
    ];
  }

  if (model === 'M/G/1') {
    const res = calculateMG1(lambda, mu, sigma);
    if (!res.isStable) {
      return [
        {
          stepNumber: 1,
          stageName: 'Phép tính 1 (P₀): Xác suất quầy rảnh rỗi',
          description: 'Hệ số sử dụng ρ = λ/μ ≥ 1. Hệ thống không ổn định.',
          latexFormula: 'P_0 = 1 - \\rho = 0',
          numericalResult: `\\rho = ${res.rho.toFixed(3)} \\ge 1`,
          explanation: 'Hệ thống bùng nổ do quá tải.',
        },
        {
          stepNumber: 2,
          stageName: 'Phép tính 2 (Lq): Chiều dài hàng chờ (Pollaczek-Khinchine)',
          description: 'Hàng chờ tiến đến vô hạn.',
          latexFormula: 'L_q \\to \\infty',
          numericalResult: 'L_q = \\infty',
          explanation: 'Hàng chờ dồn ứ vô hạn.',
        },
        {
          stepNumber: 3,
          stageName: 'Phép tính 3 (L): Số khách trong toàn hệ thống',
          description: 'Quá tải toàn diện.',
          latexFormula: 'L \\to \\infty',
          numericalResult: 'L = \\infty',
          explanation: 'Hệ thống sụp đổ.',
        },
        {
          stepNumber: 4,
          stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình',
          description: 'Thời gian chờ vô tận.',
          latexFormula: 'W_q \\to \\infty',
          numericalResult: 'W_q = \\infty',
          explanation: 'Chờ vô tận.',
        },
        {
          stepNumber: 5,
          stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình',
          description: 'Thời gian lưu lại vô tận.',
          latexFormula: 'W \\to \\infty',
          numericalResult: 'W = \\infty',
          explanation: 'Hệ thống sụp đổ.',
        },
        {
          stepNumber: 6,
          stageName: 'Phép tính 6 (Pw): Xác suất phải chờ',
          description: 'Server luôn bận 100%.',
          latexFormula: 'P_w = 100\\%',
          numericalResult: 'P_w = 100\\%',
          explanation: 'Mọi khách đến đều phải chờ.',
        },
      ];
    }

    const varS = sigma * sigma;
    return [
      {
        stepNumber: 1,
        stageName: 'Phép tính 1 (P₀): Xác suất không có khách trong hệ thống (Hệ thống rảnh)',
        description: 'Mô hình M/G/1: Lượt đến Poisson, thời gian phục vụ phân phối tổng quát G (kỳ vọng 1/μ, phương sai σ²), 1 server.',
        latexFormula: 'P_0 = 1 - \\rho = 1 - \\frac{\\lambda}{\\mu}',
        numericalResult: `P_0 = 1 - \\frac{${lambda}}{${mu}} = 1 - ${res.rho.toFixed(3)} = ${res.P0.toFixed(4)} \\quad (${(res.P0 * 100).toFixed(2)}\\%)`,
        explanation: `Quầy rảnh rỗi ${(res.P0 * 100).toFixed(2)}% thời gian. Lưu ý: Xác suất P₀ chỉ phụ thuộc vào ρ, không phụ thuộc vào phương sai σ².`,
      },
      {
        stepNumber: 2,
        stageName: 'Phép tính 2 (Lq): Số lượng khách trung bình trong hàng chờ (Công thức Pollaczek-Khinchine)',
        description: 'Công thức P-K nổi tiếng thể hiện tác động trực tiếp của độ biến thiên thời gian phục vụ (σ²) lên chiều dài hàng chờ.',
        latexFormula: 'L_q = \\frac{\\lambda^2 \\sigma^2 + \\rho^2}{2(1 - \\rho)}',
        numericalResult: `L_q = \\frac{${lambda}^2 \\times ${varS.toFixed(5)} + ${res.rho.toFixed(3)}^2}{2(1 - ${res.rho.toFixed(3)})} = ${res.Lq.toFixed(4)} \\text{ khách}`,
        explanation: `Độ biến thiên thời gian phục vụ (σ = ${sigma}) làm tăng chiều dài hàng chờ lên ${res.Lq.toFixed(2)} khách. Phương sai càng lớn thì hàng chờ càng dài!`,
      },
      {
        stepNumber: 3,
        stageName: 'Phép tính 3 (L): Số lượng khách trung bình trong toàn hệ thống',
        description: 'Tổng số khách có mặt trong hệ thống M/G/1.',
        latexFormula: 'L = L_q + \\rho = \\frac{\\lambda^2 \\sigma^2 + \\rho^2}{2(1 - \\rho)} + \\rho',
        numericalResult: `L = ${res.Lq.toFixed(4)} + ${res.rho.toFixed(4)} = ${res.L.toFixed(4)} \\text{ khách}`,
        explanation: `Trung bình có ${res.L.toFixed(2)} khách trong toàn bộ hệ thống (${res.Lq.toFixed(2)} khách đợi + ${res.rho.toFixed(2)} khách đang được phục vụ).`,
      },
      {
        stepNumber: 4,
        stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình trong hàng chờ',
        description: 'Thời gian một khách hàng phải chờ đợi trong hàng (theo Định luật Little: Wq = Lq / λ).',
        latexFormula: 'W_q = \\frac{L_q}{\\lambda} = \\frac{\\lambda \\sigma^2 + \\frac{\\rho}{\\mu}}{2(1 - \\rho)}',
        numericalResult: `W_q = \\frac{${res.Lq.toFixed(4)}}{${lambda}} = ${res.Wq.toFixed(4)} \\text{ giờ} = ${(res.Wq * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Khách hàng phải chờ trung bình ${(res.Wq * 60).toFixed(2)} phút trong hàng chờ trước khi được phục vụ.`,
      },
      {
        stepNumber: 5,
        stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình trong toàn hệ thống',
        description: 'Tổng thời gian một khách hàng lưu lại từ khi đến cho tới khi phục vụ xong.',
        latexFormula: 'W = W_q + \\frac{1}{\\mu} = \\frac{L}{\\lambda}',
        numericalResult: `W = ${res.Wq.toFixed(4)} + \\frac{1}{${mu}} = ${res.W.toFixed(4)} \\text{ giờ} = ${(res.W * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Tổng thời gian lưu lại là ${(res.W * 60).toFixed(2)} phút (gồm ${(res.Wq * 60).toFixed(2)} phút chờ và ${((1 / mu) * 60).toFixed(2)} phút phục vụ trung bình).`,
      },
      {
        stepNumber: 6,
        stageName: 'Phép tính 6 (Pw): Xác suất khách hàng khi đến phải chờ',
        description: 'Theo định lý PASTA (Poisson Arrivals See Time Averages), xác suất khách đến gặp server đang bận.',
        latexFormula: 'P_w = P(n \\ge 1) = 1 - P_0 = \\rho = \\frac{\\lambda}{\\mu}',
        numericalResult: `P_w = 1 - ${res.P0.toFixed(4)} = ${res.Pw.toFixed(4)} \\quad (${(res.Pw * 100).toFixed(2)}\\%)`,
        explanation: `Có ${(res.Pw * 100).toFixed(2)}% khả năng một khách hàng khi bước vào sẽ phải chờ vì server đang bận.`,
      },
    ];
  }

  if (model === 'M/D/1') {
    const res = calculateMD1(lambda, mu);
    if (!res.isStable) {
      return [
        {
          stepNumber: 1,
          stageName: 'Phép tính 1 (P₀): Xác suất quầy rảnh rỗi',
          description: 'Hệ số sử dụng ρ = λ/μ ≥ 1. Hệ thống không ổn định.',
          latexFormula: 'P_0 = 1 - \\rho = 0',
          numericalResult: `\\rho = ${res.rho.toFixed(3)} \\ge 1`,
          explanation: 'Hệ thống bùng nổ do quá tải.',
        },
        {
          stepNumber: 2,
          stageName: 'Phép tính 2 (Lq): Chiều dài hàng chờ M/D/1',
          description: 'Hàng chờ tiến đến vô hạn.',
          latexFormula: 'L_q \\to \\infty',
          numericalResult: 'L_q = \\infty',
          explanation: 'Hàng chờ dồn ứ vô hạn.',
        },
        {
          stepNumber: 3,
          stageName: 'Phép tính 3 (L): Số khách trong toàn hệ thống',
          description: 'Quá tải toàn diện.',
          latexFormula: 'L \\to \\infty',
          numericalResult: 'L = \\infty',
          explanation: 'Hệ thống sụp đổ.',
        },
        {
          stepNumber: 4,
          stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình',
          description: 'Thời gian chờ vô tận.',
          latexFormula: 'W_q \\to \\infty',
          numericalResult: 'W_q = \\infty',
          explanation: 'Chờ vô tận.',
        },
        {
          stepNumber: 5,
          stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình',
          description: 'Thời gian lưu lại vô tận.',
          latexFormula: 'W \\to \\infty',
          numericalResult: 'W = \\infty',
          explanation: 'Hệ thống sụp đổ.',
        },
        {
          stepNumber: 6,
          stageName: 'Phép tính 6 (Pw): Xác suất phải chờ',
          description: 'Server luôn bận 100%.',
          latexFormula: 'P_w = 100\\%',
          numericalResult: 'P_w = 100\\%',
          explanation: 'Mọi khách đến đều phải chờ do server luôn bận.',
        },
      ];
    }

    return [
      {
        stepNumber: 1,
        stageName: 'Phép tính 1 (P₀): Xác suất không có khách trong hệ thống (Hệ thống rảnh)',
        description: 'Mô hình M/D/1: Lượt đến Poisson, thời gian phục vụ cố định 1/μ (Deterministic, σ = 0), 1 server. Xác suất server hoàn toàn rảnh rỗi.',
        latexFormula: 'P_0 = 1 - \\rho = 1 - \\frac{\\lambda}{\\mu}',
        numericalResult: `P_0 = 1 - \\frac{${lambda}}{${mu}} = 1 - ${res.rho.toFixed(3)} = ${res.P0.toFixed(4)} \\quad (${(res.P0 * 100).toFixed(2)}\\%)`,
        explanation: `Thời gian phục vụ là hằng số cố định D = 1/μ = ${((1 / mu) * 60).toFixed(1)} phút. Xác suất server rảnh rỗi là ${(res.P0 * 100).toFixed(2)}%.`,
      },
      {
        stepNumber: 2,
        stageName: 'Phép tính 2 (Lq): Số lượng khách trung bình trong hàng chờ',
        description: 'Vì thời gian phục vụ cố định (σ² = 0), công thức Pollaczek-Khinchine rút gọn: Chiều dài hàng chờ Lq bằng chính xác 1/2 so với mô hình M/M/1!',
        latexFormula: 'L_q = \\frac{\\rho^2}{2(1 - \\rho)} = \\frac{\\lambda^2}{2\\mu(\\mu - \\lambda)} \\quad \\left( = \\frac{1}{2} L_{q(M/M/1)} \\right)',
        numericalResult: `L_q = \\frac{${res.rho.toFixed(3)}^2}{2(1 - ${res.rho.toFixed(3)})} = ${res.Lq.toFixed(4)} \\text{ khách}`,
        explanation: `Kết quả kinh điển: Do loại bỏ hoàn toàn sự biến thiên phục vụ (σ² = 0), hàng chờ chỉ còn ${res.Lq.toFixed(2)} khách (giảm đúng 50% so với M/M/1).`,
      },
      {
        stepNumber: 3,
        stageName: 'Phép tính 3 (L): Số lượng khách trung bình trong toàn hệ thống',
        description: 'Tổng số khách hàng trung bình có mặt trong hệ thống M/D/1.',
        latexFormula: 'L = L_q + \\rho = \\frac{\\rho^2}{2(1 - \\rho)} + \\rho = \\frac{\\rho(2 - \\rho)}{2(1 - \\rho)}',
        numericalResult: `L = ${res.Lq.toFixed(4)} + ${res.rho.toFixed(4)} = ${res.L.toFixed(4)} \\text{ khách}`,
        explanation: `Trung bình có ${res.L.toFixed(2)} khách trong toàn bộ hệ thống (${res.Lq.toFixed(2)} khách trong hàng + ${res.rho.toFixed(2)} khách đang được phục vụ).`,
      },
      {
        stepNumber: 4,
        stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình trong hàng chờ',
        description: 'Thời gian chờ trung bình trong hàng của khách (theo Định luật Little: Wq = Lq / λ).',
        latexFormula: 'W_q = \\frac{L_q}{\\lambda} = \\frac{\\rho}{2\\mu(1 - \\rho)} \\quad \\left( = \\frac{1}{2} W_{q(M/M/1)} \\right)',
        numericalResult: `W_q = \\frac{${res.Lq.toFixed(4)}}{${lambda}} = ${res.Wq.toFixed(4)} \\text{ giờ} = ${(res.Wq * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Khách hàng chỉ phải chờ trung bình ${(res.Wq * 60).toFixed(2)} phút trong hàng (bằng 1/2 so với mô hình M/M/1 cùng thông số).`,
      },
      {
        stepNumber: 5,
        stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình trong toàn hệ thống',
        description: 'Tổng thời gian một khách hàng lưu lại từ khi đến cho tới khi kết thúc phục vụ rời đi.',
        latexFormula: 'W = W_q + \\frac{1}{\\mu} = \\frac{L}{\\lambda}',
        numericalResult: `W = ${res.Wq.toFixed(4)} + \\frac{1}{${mu}} = ${res.W.toFixed(4)} \\text{ giờ} = ${(res.W * 60).toFixed(2)} \\text{ phút}`,
        explanation: `Tổng thời gian lưu lại trung bình là ${(res.W * 60).toFixed(2)} phút (gồm ${(res.Wq * 60).toFixed(2)} phút chờ hàng và đúng ${((1 / mu) * 60).toFixed(2)} phút phục vụ cố định).`,
      },
      {
        stepNumber: 6,
        stageName: 'Phép tính 6 (Pw): Xác suất khách hàng khi đến phải chờ',
        description: 'Xác suất một khách hàng khi đến thấy server đang bận phục vụ và buộc phải vào hàng chờ.',
        latexFormula: 'P_w = P(n \\ge 1) = 1 - P_0 = \\rho = \\frac{\\lambda}{\\mu}',
        numericalResult: `P_w = 1 - ${res.P0.toFixed(4)} = ${res.Pw.toFixed(4)} \\quad (${(res.Pw * 100).toFixed(2)}\\%)`,
        explanation: `Dù thời gian phục vụ là hằng số, nhưng lượt đến Poisson ngẫu nhiên (PASTA) nên xác suất phải chờ vẫn bằng tỷ lệ bận server (${(res.Pw * 100).toFixed(2)}%).`,
      },
    ];
  }

  // Model: M/G/k không có hàng chờ (M/G/k/k Loss System / Erlang-B)
  const res = calculateMGkLoss(lambda, mu, k);
  return [
    {
      stepNumber: 1,
      stageName: 'Phép tính 1 (P₀): Xác suất không có khách trong hệ thống (Tất cả quầy đều rảnh)',
      description: `Mô hình M/G/k không có hàng chờ (M/G/k/k Loss System): Lượt đến Poisson λ = ${lambda} khách/giờ, thời gian phục vụ tổng quát G (trung bình 1/μ = ${((1 / mu) * 60).toFixed(1)} phút), k = ${k} servers, dung lượng hàng chờ = 0. Tải cung cấp a = λ/μ.`,
      latexFormula: 'P_0 = \\frac{1}{\\sum_{n=0}^k \\frac{a^n}{n!}} \\quad \\text{với } a = \\frac{\\lambda}{\\mu}',
      numericalResult: `a = \\frac{${lambda}}{${mu}} = ${res.a.toFixed(3)} \\text{ Erlang} \\implies P_0 = ${res.P0.toFixed(4)} \\quad (${(res.P0 * 100).toFixed(2)}\\%)`,
      explanation: `Đặc tính bất biến Erlang (Insensitivity Property): Xác suất P₀ chỉ phụ thuộc vào tỷ lệ tải a = λ/μ = ${res.a.toFixed(2)}, hoàn toàn độc lập với hình dạng phân phối hay phương sai σ² của thời gian phục vụ G.`,
    },
    {
      stepNumber: 2,
      stageName: 'Phép tính 2 (Lq): Số lượng khách trung bình trong hàng chờ',
      description: 'Vì hệ thống được thiết kế không có hàng chờ (sức chứa hàng chờ bằng 0, khách đến khi k quầy bận sẽ bị từ chối / rời đi ngay).',
      latexFormula: 'L_q = 0 \\text{ (Không có vùng đệm chờ)}',
      numericalResult: 'L_q = 0 \\text{ khách}',
      explanation: 'Do hệ thống không lưu giữ khách xếp hàng, chiều dài hàng chờ Lq luôn bằng 0 tuyệt đối.',
    },
    {
      stepNumber: 3,
      stageName: 'Phép tính 3 (L): Số lượng khách trung bình trong toàn hệ thống (Số quầy đang bận)',
      description: 'Tổng số khách trung bình có mặt trong hệ thống bằng chính số quầy đang bận phục vụ tại trạng thái dừng.',
      latexFormula: 'L = a(1 - P_k) = \\frac{\\lambda(1 - P_k)}{\\mu} = \\frac{\\lambda_{eff}}{\\mu}',
      numericalResult: `L = ${res.a.toFixed(3)} \\times (1 - ${res.Pk.toFixed(4)}) = ${res.L.toFixed(4)} \\text{ khách}`,
      explanation: `Thông lượng khách thực tế được tiếp nhận là λ_eff = λ(1 - Pk) = ${res.lambdaEff.toFixed(2)} khách/giờ. Trung bình có ${res.L.toFixed(2)} / ${k} quầy đang bận phục vụ.`,
    },
    {
      stepNumber: 4,
      stageName: 'Phép tính 4 (Wq): Thời gian chờ trung bình trong hàng chờ',
      description: 'Thời gian chờ trong hàng của khách hàng trước khi vào quầy phục vụ.',
      latexFormula: 'W_q = 0 \\text{ (Không có hàng chờ)}',
      numericalResult: 'W_q = 0 \\text{ giờ} = 0 \\text{ phút}',
      explanation: 'Khách hàng đến khi có quầy rảnh sẽ được phục vụ ngay lập tức, không mất thời gian đợi.',
    },
    {
      stepNumber: 5,
      stageName: 'Phép tính 5 (W): Thời gian lưu lại trung bình trong toàn hệ thống',
      description: 'Thời gian trung bình của một khách hàng được tiếp nhận từ khi vào đến khi rời hệ thống.',
      latexFormula: 'W = W_q + \\frac{1}{\\mu} = 0 + \\frac{1}{\\mu} = \\frac{1}{\\mu} = \\frac{L}{\\lambda_{eff}}',
      numericalResult: `W = \\frac{1}{${mu}} = ${res.W.toFixed(4)} \\text{ giờ} = ${(res.W * 60).toFixed(2)} \\text{ phút}`,
      explanation: `Thời gian lưu lại đúng bằng thời gian phục vụ trung bình tại quầy: ${(res.W * 60).toFixed(2)} phút (thỏa mãn chính xác Định luật Little: W = L / λ_eff).`,
    },
    {
      stepNumber: 6,
      stageName: 'Phép tính 6 (Pw hay P_loss): Xác suất mất khách (Công thức Erlang B)',
      description: `Xác suất một khách hàng khi đến gặp toàn bộ k = ${k} quầy đều đang bận và bị từ chối phục vụ (Blocking probability / Erlang-B).`,
      latexFormula: 'P_w = P_k = B(k, a) = \\frac{\\frac{a^k}{k!}}{\\sum_{n=0}^k \\frac{a^n}{n!}}',
      numericalResult: `P_w = P_{${k}} = B(${k}, ${res.a.toFixed(2)}) = ${res.Pk.toFixed(4)} \\quad (${(res.Pk * 100).toFixed(2)}\\%)`,
      explanation: `Xác suất bị chặn / mất khách là ${(res.Pk * 100).toFixed(2)}%. Khách hàng bị mất trung bình mỗi giờ là: λ × Pk = ${(lambda * res.Pk).toFixed(2)} khách/giờ.`,
    },
  ];
}

export function calculateOptimalServers(
  lambda: number,
  mu: number,
  costService: number,
  costWaiting: number,
  maxServers: number = 8
) {
  const minServers = Math.max(1, Math.ceil(lambda / mu));
  const results = [];

  for (let c = minServers; c <= maxServers; c++) {
    const metrics = calculateMMC(lambda, mu, c);
    if (!metrics.isStable) continue;

    const costS = costService * c;
    const costW = costWaiting * metrics.L;
    const totalCost = costS + costW;

    results.push({
      servers: c,
      utilization: metrics.rho,
      L: metrics.L,
      Lq: metrics.Lq,
      W: metrics.W,
      Wq: metrics.Wq,
      costService: costS,
      costWaiting: costW,
      totalCost,
      isOptimal: false,
    });
  }

  if (results.length > 0) {
    let minCost = Infinity;
    let minIdx = 0;
    results.forEach((r, idx) => {
      if (r.totalCost < minCost) {
        minCost = r.totalCost;
        minIdx = idx;
      }
    });
    results[minIdx].isOptimal = true;
  }

  return results;
}
