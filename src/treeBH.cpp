#include <Rcpp.h>
#include <string>
#include <map>
#include <set>
#include <algorithm>
using namespace Rcpp;
using namespace std;

//' Check if groups are properly nested in hierarchy
// [[Rcpp::export]]
std::string checkNested(const CharacterMatrix& groups, int L) {
  int nrow = groups.nrow();
  int ncol = groups.ncol();

  // Validate L to avoid out-of-bounds column access
  if (L < 1 || L > ncol) {
    return "Invalid hierarchy level L: L must be between 1 and the number of columns in 'groups'.";
  }
  // Step 1: Check hierarchy nesting
  for (int cur_level = 1; cur_level < L; cur_level++) {

    // Check if any child has multiple parents
    map<string, set<string>> child_to_parents;
    for (int i = 0; i < nrow; i++) {
      string child = as<string>(groups(i, cur_level));
      string parent = as<string>(groups(i, cur_level - 1));

      child_to_parents[child].insert(parent);

      // If child has more than one parent, return error
      if (child_to_parents[child].size() > 1) {
        return "Groups must be nested within hierarchy";
      }
    }
  }

  // Step 2: Check if lowest level corresponds to individual hypotheses
  set<string> lowest_level_values;
  for (int i = 0; i < nrow; i++) {
    lowest_level_values.insert(as<string>(groups(i, L - 1)));
  }

  if (lowest_level_values.size() != nrow) {
    return "Assumption is that lowest level in tree corresponds to individual hypotheses";
  }

  // All checks passed
  return "";
}

// Calculate Simes p-value
double simesPvalue(const NumericVector& pvals) {
  // Remove NA values
  NumericVector clean_pvals;
  for (int i = 0; i < pvals.size(); i++) {
    if (!NumericVector::is_na(pvals[i])) {
      clean_pvals.push_back(pvals[i]);
    }
  }

  // If no valid p-values, return NA
  if (clean_pvals.size() == 0) {
    return NA_REAL;
  }

  // Sort p-values
  std::sort(clean_pvals.begin(), clean_pvals.end());

  // Calculate Simes p-value
  double min_val = R_PosInf;
  int n = clean_pvals.size();

  for (int i = 0; i < n; i++) {
    double adjusted = clean_pvals[i] / (i + 1);
    if (adjusted < min_val) {
      min_val = adjusted;
    }
  }

  double pvalue = n * min_val;
  if (pvalue < 0.0) {
    pvalue = 0.0;
  } else if (pvalue > 1.0) {
    pvalue = 1.0;
  }
  return pvalue;
}

// Calculate Fisher's pvalue
double fisherPvalue(const NumericVector& pvals) {
  // Remove NA values
  NumericVector clean_pvals;
  for (int i = 0; i < pvals.size(); i++) {
    if (!NumericVector::is_na(pvals[i])) {
      clean_pvals.push_back(pvals[i]);
    }
  }
  
  // If no valid p-values, return NA
  if (clean_pvals.size() == 0) {
    return NA_REAL;
  }
  
  // Calculate Fisher's statistic: -2 * sum(log(p))
  double sum_log = 0.0;
  for (int i = 0; i < clean_pvals.size(); i++) {
    sum_log += log(clean_pvals[i]);
  }
  double statistic = -2 * sum_log;
  
  // Return p-value from chi-squared distribution with 2k degrees of freedom
  // Using R's built-in function
  return R::pchisq(statistic, 2 * clean_pvals.size(), false, false);
}

//' Hierarchical Group P-value Calculation
//' The getGroupPvalues function computes aggregated p-values for each group in a hierarchical 
//' structure by combining p-values from lower levels using specified aggregation methods.
//' @param pvals NumericVector containing p-values for individual hypotheses
//' @param groups CharacterMatrix defining the group hierarchy structure
//' @param N Integer specifying the number of hypotheses
//' @param L Integer specifying the number of levels in the hierarchy
//' @param test CharacterVector specifying the aggregation method(s)
// [[Rcpp::export]]
NumericMatrix getGroupPvalues(const NumericVector& pvals, const CharacterMatrix& groups,
                              int N, int L, CharacterVector test) {
  // Initialize group_pvals matrix
  NumericMatrix group_pvals(N, L);
  
  // Fill in initial values
  std::fill(group_pvals.begin(), group_pvals.end(), NA_REAL);
  
  // Set initial values for group_pvals at the lowest level
  for (int i = 0; i < N; i++) {
    group_pvals(i, L-1) = pvals[i];
  }
  
  // Ensure test vector is of correct length
  CharacterVector test_methods = test;
  if (test_methods.size() == 1) {
    test_methods = CharacterVector(L-1, test[0]);
  }
  
  // Pre-compute group indices for each level
  std::vector<std::map<string, std::vector<int>>> level_groups(L);
  
  // Build group index maps for each level
  for (int level = 0; level < L; level++) {
    for (int i = 0; i < N; i++) {
      string group = as<string>(groups(i, level));
      level_groups[level][group].push_back(i);
    }
  }
  
  // Step 3: Aggregate group-level p-values
  for (int cur_level = L - 2; cur_level >= 0; cur_level--) {
    // Get unique groups at the current level
    const auto& groups_at_level = level_groups[cur_level];
    
    // Cache test method for this level
    string method = as<string>(test_methods[cur_level]);
    
    // Process each group
    for (const auto& group_entry : groups_at_level) {
      const std::vector<int>& cur_inds = group_entry.second;
      
      // Skip empty groups (shouldn't happen but check anyway)
      if (cur_inds.empty()) continue;
      
      // Get child p-values
      NumericVector child_pvals(cur_inds.size());
      for (size_t i = 0; i < cur_inds.size(); i++) {
        child_pvals[i] = group_pvals(cur_inds[i], cur_level + 1);
      }
      
      // Calculate aggregated p-value
      double agg_p;
      if (method == "simes") {
        agg_p = simesPvalue(child_pvals);
      } else if (method == "fisher") {
        agg_p = fisherPvalue(child_pvals);
      } else {
        stop("Options for parameter test are 'simes' and 'fisher'");
      }
      
      // Set the aggregated p-value ONLY for the first index in the group
      group_pvals(cur_inds[0], cur_level) = agg_p;
    }
  }
  
  // Return only the group_pvals matrix
  return group_pvals;
}

// Function to compute qvalues (simplified version using BH method)
NumericVector computeQValues(const NumericVector& pvals) {
  int n = pvals.size();
  if (n == 0) return NumericVector(0);
  
  // Create index vector and sort by p-value
  IntegerVector idx = seq(0, n-1);
  std::sort(idx.begin(), idx.end(), [&pvals](int i, int j) {
    return pvals[i] < pvals[j];
  });
  
  // Calculate adjusted p-values
  NumericVector adjusted(n);
  adjusted[idx[n-1]] = std::min(1.0, pvals[idx[n-1]]);
  
  for (int i = n-2; i >= 0; i--) {
    adjusted[idx[i]] = std::min(adjusted[idx[i+1]], pvals[idx[i]] * n / (i + 1));
  }
  
  return adjusted;
}

//' Hierarchical Group Selection for Multiple Testing
//' The getGroupSelections function applies hierarchical multiple testing correction.
//' @param group_pvals NumericMatrix containing aggregated p-values
//' @param groups CharacterMatrix defining the group hierarchy structure
//' @param q NumericVector containing the q-value thresholds for each level
//' @param N Integer specifying the number of hypotheses
//' @param L Integer specifying the number of levels in the hierarchy
// [[Rcpp::export]]
NumericMatrix getGroupSelections(const NumericMatrix& group_pvals, 
                               const CharacterMatrix& groups,
                               const NumericVector& q, int N, int L) {
  // Validate input dimensions to prevent out-of-bounds access
  if (L <= 0) {
    Rcpp::stop("getGroupSelections: L must be at least 1.");
  }
  if (group_pvals.nrow() != N || group_pvals.ncol() != L) {
    Rcpp::stop("getGroupSelections: N and L must match group_pvals dimensions (nrow, ncol).");
  }
  if (groups.nrow() != N || groups.ncol() != L) {
    Rcpp::stop("getGroupSelections: N and L must match groups dimensions (nrow, ncol).");
  }
  if (q.size() != L) {
    Rcpp::stop("getGroupSelections: length(q) must be equal to L.");
  }
  NumericMatrix sel(N, L);
  NumericMatrix q_adj(N, L);
  
  // Initialize selection matrix to 0 (no selections initially)
  std::fill(sel.begin(), sel.end(), 0.0);
  
  // Fill q_adj with NA
  std::fill(q_adj.begin(), q_adj.end(), NA_REAL);
  
  // Set initial value for q_adj
  for (int i = 0; i < N; i++) {
    q_adj(i, 0) = q[0];
  }
  
  // Iterate through each level
  for (int cur_level = 0; cur_level < L; cur_level++) {
    // Check if we should break
    if (cur_level > 0) {
      bool any_selected = false;
      for (int i = 0; i < N; i++) {
        if (sel(i, cur_level - 1) == 1) {
          any_selected = true;
          break;
        }
      }
      if (!any_selected) break;
    }
    
    // Get selected parents from previous level
    std::vector<int> sel_prev;
    if (cur_level == 0) {
      // First level: select all
      sel_prev.push_back(0); // We'll use a dummy index 0 for the root
    } else {
      // Later levels: get indices of selected items from previous level
      for (int i = 0; i < N; i++) {
        if (sel(i, cur_level - 1) == 1) {
          sel_prev.push_back(i);
        }
      }
    }
    
    // Process each selected parent
    for (int parent_sel : sel_prev) {
      std::vector<int> child_inds;
      int parent_group_ind;
      
      if (cur_level == 0) {
        // First level: include all non-NA indices
        parent_group_ind = 0;
        for (int i = 0; i < N; i++) {
          if (!NumericVector::is_na(group_pvals(i, cur_level))) {
            child_inds.push_back(i);
          }
        }
      } else {
        // Later levels: include children of the current parent group
        string parent_group_num = as<string>(groups(parent_sel, cur_level - 1));
        
        // Find the first index with this parent group
        for (int i = 0; i < N; i++) {
          if (as<string>(groups(i, cur_level - 1)) == parent_group_num) {
            parent_group_ind = i;
            break;
          }
        }
        
        // Find all indices with this parent group and non-NA p-values
        for (int i = 0; i < N; i++) {
          if (as<string>(groups(i, cur_level - 1)) == parent_group_num && 
              !NumericVector::is_na(group_pvals(i, cur_level))) {
            child_inds.push_back(i);
          }
        }
      }
      
      // Skip if no children found
      if (child_inds.empty()) continue;
      
      // Select children based on p-values and thresholds
      std::vector<int> sel_ind_within_group;
      
      if (child_inds.size() > 1) {
        // Multiple children: use q-value method
        NumericVector child_pvals(child_inds.size());
        for (size_t i = 0; i < child_inds.size(); i++) {
          child_pvals[i] = group_pvals(child_inds[i], cur_level);
        }
        
        // Compute q-values
        NumericVector qvals = computeQValues(child_pvals);
        
        // Select children with q-values <= threshold
        for (size_t i = 0; i < child_inds.size(); i++) {
          if (qvals[i] <= q_adj(parent_group_ind, cur_level)) {
            sel_ind_within_group.push_back(i);
          }
        }
      } else {
        // Single child: use direct p-value comparison
        if (group_pvals(child_inds[0], cur_level) <= q_adj(parent_group_ind, cur_level)) {
          sel_ind_within_group.push_back(0);
        }
      }
      
      // Update selection matrix
      for (int idx : sel_ind_within_group) {
        sel(child_inds[idx], cur_level) = 1;
      }
      
      // Update q_adj for the next level
      if (cur_level < L - 1) {
        double R_parent_sel = sel_ind_within_group.size();
        double n_parent_sel = child_inds.size();
        
        if (R_parent_sel > 0) {
          for (int idx : sel_ind_within_group) {
            q_adj(child_inds[idx], cur_level + 1) = q[cur_level + 1] * 
              q_adj(parent_sel, cur_level) / q[cur_level] * 
              R_parent_sel / n_parent_sel;
          }
        }
      }
    }
  }
  
  return sel;
}
