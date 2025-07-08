# A Quadratic Programming Algorithm with $O(n^3)$ Time Complexity 
This is the code repository for https://arxiv.org/pdf/2507.04515

(1) demonstrates how to transform a general strictly convex QP, a Lasso problem, and a support vector machine (or regression) into a Box-QP, highlighting the broad applicability of the approach.

(2) develops a unified and simple proof framework for feasible IPM algorithms with \textit{exact Newton step} and \textit{approximated Newton step}, respectively, resulting in $O(n^{3.5})$ and $O(n^3)$ time complexity.

(3)  proves that the proposed $O(n^3)$-time-complexity algorithm has an \textit{exact} and \textit{data-independent} number of iterations   
$\mathcal{N}_\mathrm{iter}= \!\left\lceil\frac{\log\!\left(\frac{2n+\alpha\sqrt{2n}}{\epsilon}\right)}{-\log\!\left(1-\frac{\beta}{\sqrt{2n}}\right)}\right\rceil$
and the \textit{data-independent} total number of rank-1 updates bounded by $\mathcal{N}_{\mathrm{rank}-1}\leq \! \left\lceil\frac{4\eta(\mathcal{N}_\mathrm{iter}-1)\sqrt{n}}{\left(1-\eta\right)\log(1+\delta)}\right\rceil.$
where $\alpha,\beta,\eta,\delta$ are \textit{data-independent} constants, thereby being able to offer an execution time certificate for real-time optimization-based applications;

(4) shows that the proposed $O(n^3)$-time-complexity algorithm is simple to implement and matrix-free in the iterative procedures.

