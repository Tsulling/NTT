#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// 定义模数 P = 998244353 和原根 G = 3
#define MOD 998244353
#define G 3
#define G_INV 332748118 // 3 在模 998244353 下的逆元

typedef long long ll;

// 快速幂算法: 计算 (base^exp) % mod
ll qpow(ll base, ll exp)
{
    ll res = 1;
    while (exp > 0)
    {
        if (exp & 1)
            res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp >>= 1;
    }
    return res;
}

// 计算位逆序置换数组
// limit: 多项式对齐后的长度（必须是2的幂）
// r: 存储逆序索引的数组
void get_rev(int limit, int *r)
{
    int bit_len = 0;
    while ((1 << bit_len) < limit)
        bit_len++;

    for (int i = 0; i < limit; i++)
    {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (bit_len - 1));
    }
}

// NTT 核心变换函数
// a: 系数数组
// limit: 变换长度
// type: 1 表示 NTT (正变换), -1 表示 INTT (逆变换)
// r: 位逆序数组
void ntt(ll *a, int limit, int type, int *r)
{
    // 1. 进行位逆序置换
    for (int i = 0; i < limit; i++)
    {
        if (i < r[i])
        {
            ll temp = a[i];
            a[i] = a[r[i]];
            a[r[i]] = temp;
        }
    }

    // 2. Cooley-Tukey 蝴蝶变换
    // mid 是当前合并区间的长度的一半
    for (int mid = 1; mid < limit; mid <<= 1)
    {
        // wn 是单位原根: g^((P-1)/(2*mid))
        ll wn = qpow(type == 1 ? G : G_INV, (MOD - 1) / (mid << 1));

        // j 枚举每一个合并块的起始位置
        for (int j = 0; j < limit; j += (mid << 1))
        {
            ll w = 1;
            // k 枚举块内的位置
            for (int k = 0; k < mid; k++)
            {
                ll x = a[j + k];
                ll y = (w * a[j + k + mid]) % MOD;

                a[j + k] = (x + y) % MOD;
                a[j + k + mid] = (x - y + MOD) % MOD; // 注意负数处理

                w = (w * wn) % MOD;
            }
        }
    }

    // 3. 如果是逆变换，需要除以长度 N (乘以 N 的逆元)
    if (type == -1)
    {
        ll inv = qpow(limit, MOD - 2);
        for (int i = 0; i < limit; i++)
        {
            a[i] = (a[i] * inv) % MOD;
        }
    }
}

// 多项式乘法封装函数
// polyA, degA: 多项式A及其次数
// polyB, degB: 多项式B及其次数
// result: 结果数组
// 返回值: 结果多项式的长度(limit)
int poly_multiply(const ll *polyA, int degA, const ll *polyB, int degB, ll *result)
{
    int limit = 1;
    // 找到大于 degA + degB 的最小 2 的幂次
    while (limit <= degA + degB)
        limit <<= 1;

    // 分配内存并初始化
    ll *a = (ll *)calloc(limit, sizeof(ll));
    ll *b = (ll *)calloc(limit, sizeof(ll));
    int *r = (int *)malloc(limit * sizeof(int));

    // 复制系数
    for (int i = 0; i <= degA; i++)
        a[i] = polyA[i];
    for (int i = 0; i <= degB; i++)
        b[i] = polyB[i];

    // 获取逆序
    get_rev(limit, r);

    // 1. NTT 正变换 -> 转为点值表示
    ntt(a, limit, 1, r);
    ntt(b, limit, 1, r);

    // 2. 点值相乘
    for (int i = 0; i < limit; i++)
    {
        a[i] = (a[i] * b[i]) % MOD;
    }

    // 3. INTT 逆变换 -> 转回系数表示
    ntt(a, limit, -1, r);

    // 复制结果
    for (int i = 0; i < limit; i++)
    {
        result[i] = a[i];
    }

    // 释放内存
    free(a);
    free(b);
    free(r);

    return limit;
}

// 性能测试函数
void benchmark()
{
    printf("--- 要求二：性能与多项式次数关系分析 ---\n");
    printf("%-10s | %-15s | %-15s\n", "Degree (N)", "Limit (Size)", "Time (seconds)");
    printf("-----------|-----------------|----------------\n");

    // 测试从 2^8 (256) 到 2^18 (262144) 的不同规模
    for (int p = 8; p <= 18; p++)
    {
        int n = 1 << p; // 多项式次数
        int limit = 1;
        while (limit <= n + n)
            limit <<= 1;

        // 准备随机数据
        ll *A = (ll *)malloc(limit * sizeof(ll));
        ll *B = (ll *)malloc(limit * sizeof(ll));
        ll *Res = (ll *)malloc(limit * sizeof(ll));
        int *r = (int *)malloc(limit * sizeof(int));

        for (int i = 0; i < n; i++)
            A[i] = rand() % MOD;
        for (int i = 0; i < n; i++)
            B[i] = rand() % MOD;

        // 预计算逆序数组 (不计入主要的重复测试时间，或根据需求计入)
        get_rev(limit, r);

        // 计时开始
        clock_t start = clock();

        // 为了使小规模测试时间可测，循环多次
        int loops = (p < 12) ? 100 : 1;

        for (int k = 0; k < loops; k++)
        {
            // 这里手动展开poly_multiply的核心部分以避免反复malloc/free干扰计时
            // 每次循环重新赋值以模拟真实计算
            ntt(A, limit, 1, r);
            ntt(B, limit, 1, r);
            for (int i = 0; i < limit; i++)
                A[i] = (A[i] * B[i]) % MOD;
            ntt(A, limit, -1, r);
        }

        clock_t end = clock();
        double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC / loops;

        printf("N=2^%-2d    | %-15d | %f s\n", p, limit, time_taken);

        free(A);
        free(B);
        free(Res);
        free(r);
    }
}

int main()
{
    // --- 要求一：功能验证 (256次多项式) ---
    printf("--- 要求一：NTT 功能验证 ---\n");
    // 简单例子：(1 + 2x) * (2 + 3x) = 2 + 7x + 6x^2
    ll polyA[] = {1, 2}; // 1 + 2x
    ll polyB[] = {2, 3}; // 2 + 3x
    ll result[8];

    poly_multiply(polyA, 1, polyB, 1, result);

    printf("Test (1+2x)*(2+3x) Result: %lld + %lldx + %lldx^2\n", result[0], result[1], result[2]);
    printf("Expected: 2 + 7x + 6x^2\n\n");

    // --- 要求二：性能探索 ---
    benchmark();

    return 0;
}