import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import re
from itertools import combinations, chain

class NormalizationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Database Normalization Tool")
        self.root.geometry("900x700")
        self.root.configure(bg="#f0f0f0")
        
        # Create style
        style = ttk.Style()
        style.configure("TFrame", background="#f0f0f0")
        style.configure("TLabel", background="#f0f0f0", font=("Arial", 10))
        style.configure("TButton", font=("Arial", 10))
        style.configure("Header.TLabel", font=("Arial", 12, "bold"))
        
        # Main frame
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Input section
        input_frame = ttk.Frame(main_frame)
        input_frame.pack(fill=tk.X, pady=10)
        
        # Schema input
        ttk.Label(input_frame, text="Relational Schema (comma separated attributes):", style="Header.TLabel").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.schema_entry = tk.Entry(input_frame, width=70, font=("Arial", 10))
        self.schema_entry.grid(row=1, column=0, sticky=tk.W+tk.E, pady=2)
        self.schema_entry.insert(0, "A, B, C, D, E")
        
        # Functional dependencies input
        ttk.Label(input_frame, text="Functional Dependencies (one per line, format: A -> B or A,B -> C,D):", style="Header.TLabel").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.fd_text = scrolledtext.ScrolledText(input_frame, width=70, height=5, font=("Arial", 10))
        self.fd_text.grid(row=3, column=0, sticky=tk.W+tk.E, pady=2)
        self.fd_text.insert(tk.END, "A -> B\nB -> C\nD -> E\nA,D -> C")
        
        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=5)
        
        # Analysis buttons
        ttk.Button(button_frame, text="Find Attribute Closures", command=self.find_closures).grid(row=0, column=0, padx=5)
        ttk.Button(button_frame, text="Find Candidate Keys", command=self.find_candidate_keys).grid(row=0, column=1, padx=5)
        ttk.Button(button_frame, text="Normalize Schema", command=self.normalize_schema).grid(row=0, column=2, padx=5)
        ttk.Button(button_frame, text="Clear Results", command=self.clear_results).grid(row=0, column=3, padx=5)
        
        # Output section
        output_frame = ttk.Frame(main_frame)
        output_frame.pack(fill=tk.BOTH, expand=True, pady=10)
        
        ttk.Label(output_frame, text="Results:", style="Header.TLabel").pack(anchor=tk.W)
        self.output_text = scrolledtext.ScrolledText(output_frame, width=80, height=25, font=("Courier New", 10))
        self.output_text.pack(fill=tk.BOTH, expand=True)
        
        # Initialize state variables
        self.attributes = []
        self.fds = []
        self.candidate_keys = []
        
    def parse_input(self):
        """Parse user input for attributes and functional dependencies."""
        try:
            # Parse schema
            schema_text = self.schema_entry.get().strip()
            self.attributes = [attr.strip() for attr in schema_text.split(',')]
            
            # Filter out empty strings
            self.attributes = [attr for attr in self.attributes if attr]
            
            if not self.attributes:
                raise ValueError("No attributes found in schema")
                
            # Remove duplicates while maintaining order
            self.attributes = list(dict.fromkeys(self.attributes))
            
            # Parse functional dependencies
            fd_text = self.fd_text.get("1.0", tk.END).strip()
            fd_lines = fd_text.split('\n')
            
            self.fds = []
            for line in fd_lines:
                if not line.strip():
                    continue
                    
                if "->" not in line:
                    raise ValueError(f"Invalid functional dependency format: {line}")
                    
                lhs, rhs = [side.strip() for side in line.split('->')]
                
                lhs_attrs = [attr.strip() for attr in lhs.split(',')]
                rhs_attrs = [attr.strip() for attr in rhs.split(',')]
                
                # Filter empty strings
                lhs_attrs = [attr for attr in lhs_attrs if attr]
                rhs_attrs = [attr for attr in rhs_attrs if attr]
                
                if not lhs_attrs or not rhs_attrs:
                    raise ValueError(f"Empty left or right side in FD: {line}")
                    
                # Check if attributes are in schema
                for attr in lhs_attrs + rhs_attrs:
                    if attr not in self.attributes:
                        raise ValueError(f"Attribute '{attr}' in FD is not in schema")
                
                self.fds.append((frozenset(lhs_attrs), frozenset(rhs_attrs)))
                
            return True
            
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return False
    
    def compute_closure(self, attrs):
        """Compute the closure of a set of attributes."""
        closure = set(attrs)
        changed = True
        
        while changed:
            changed = False
            for lhs, rhs in self.fds:
                if lhs.issubset(closure) and not rhs.issubset(closure):
                    closure.update(rhs)
                    changed = True
                    
        return closure
    
    def find_closures(self):
        """Find and display the closure of each attribute and attribute subset."""
        if not self.parse_input():
            return
            
        self.output_text.delete("1.0", tk.END)
        self.output_text.insert(tk.END, "=== ATTRIBUTE CLOSURES ===\n\n")
        
        # For each individual attribute
        for attr in self.attributes:
            closure = self.compute_closure([attr])
            self.output_text.insert(tk.END, f"{attr}⁺ = {{{', '.join(sorted(closure))}}}\n")
            
        self.output_text.insert(tk.END, "\n=== ATTRIBUTE SET CLOSURES ===\n\n")
        
        # For pairs of attributes
        for i, attr1 in enumerate(self.attributes):
            for attr2 in self.attributes[i+1:]:
                attrs = [attr1, attr2]
                closure = self.compute_closure(attrs)
                self.output_text.insert(tk.END, f"{{{', '.join(attrs)}}}⁺ = {{{', '.join(sorted(closure))}}}\n")
    
    def find_candidate_keys(self):
        """Find and display all candidate keys for the schema."""
        if not self.parse_input():
            return
            
        self.output_text.delete("1.0", tk.END)
        self.output_text.insert(tk.END, "=== FINDING CANDIDATE KEYS ===\n\n")
        
        self.output_text.insert(tk.END, "Step 1: Identify all attributes that never appear on the right-hand side of any FD\n")
        
        # Find attributes that don't appear on RHS of any FD
        rhs_attrs = set()
        for _, rhs in self.fds:
            rhs_attrs.update(rhs)
            
        lhs_only_attrs = [attr for attr in self.attributes if attr not in rhs_attrs]
        
        if lhs_only_attrs:
            self.output_text.insert(tk.END, f"Attributes that don't appear on any RHS: {', '.join(lhs_only_attrs)}\n")
            self.output_text.insert(tk.END, "These must be part of any candidate key.\n\n")
        else:
            self.output_text.insert(tk.END, "No attributes found that appear only on LHS.\n\n")
        
        # Try to find minimal sets that determine all attributes
        self.output_text.insert(tk.END, "Step 2: Find minimal sets of attributes that determine all other attributes\n\n")
        
        self.candidate_keys = []
        
        # Start with the known LHS-only attributes
        base_set = set(lhs_only_attrs)
        remaining_attrs = [attr for attr in self.attributes if attr not in base_set]
        
        if base_set:
            self.output_text.insert(tk.END, f"Starting with required attributes: {', '.join(sorted(base_set))}\n")
        else:
            self.output_text.insert(tk.END, "No required attributes for candidate keys.\n")
        
        # If base set already determines all attributes, it's a candidate key
        if base_set and self.compute_closure(base_set) == set(self.attributes):
            self.candidate_keys.append(frozenset(base_set))
            self.output_text.insert(tk.END, f"Found candidate key: {{{', '.join(sorted(base_set))}}}\n")
        else:
            # Try adding combinations of remaining attributes until we find candidate keys
            self.output_text.insert(tk.END, "Testing combinations of attributes to find candidate keys...\n")
            
            # Generate powersets of remaining attributes
            for r in range(1, len(remaining_attrs) + 1):
                for combo in combinations(remaining_attrs, r):
                    current_set = base_set.union(combo)
                    closure = self.compute_closure(current_set)
                    
                    if closure == set(self.attributes):
                        # Check if it's minimal
                        is_minimal = True
                        for existing_key in self.candidate_keys:
                            if existing_key.issubset(current_set):
                                is_minimal = False
                                break
                                
                        if is_minimal:
                            # Check if any attribute can be removed
                            for attr in list(combo):
                                test_set = current_set - {attr}
                                if self.compute_closure(test_set) == set(self.attributes):
                                    is_minimal = False
                                    break
                            
                            if is_minimal:
                                self.candidate_keys.append(frozenset(current_set))
                                self.output_text.insert(tk.END, f"Found candidate key: {{{', '.join(sorted(current_set))}}}\n")
        
        # Display final results
        if self.candidate_keys:
            self.output_text.insert(tk.END, "\n=== CANDIDATE KEYS ===\n")
            for i, key in enumerate(self.candidate_keys, 1):
                self.output_text.insert(tk.END, f"Key {i}: {{{', '.join(sorted(key))}}}\n")
        else:
            self.output_text.insert(tk.END, "\nNo candidate keys found!\n")
    
    def normalize_schema(self):
        """Perform normalization from 1NF to BCNF with steps."""
        if not self.parse_input():
            return
            
        # First find candidate keys if not already done
        if not self.candidate_keys:
            self.find_candidate_keys()
            
        self.output_text.delete("1.0", tk.END)
        self.output_text.insert(tk.END, "=== NORMALIZATION PROCESS ===\n\n")
        
        # 1NF
        self.output_text.insert(tk.END, "STEP 1: FIRST NORMAL FORM (1NF)\n")
        self.output_text.insert(tk.END, "A relation is in 1NF if it contains no repeating groups or arrays.\n")
        self.output_text.insert(tk.END, "Assuming the input schema is already in 1NF since we're working with atomic attributes.\n\n")
        
        # Original schema in 1NF
        schema_1nf = [set(self.attributes)]
        self.output_text.insert(tk.END, f"1NF Schema: R({', '.join(self.attributes)})\n")
        self.output_text.insert(tk.END, f"Functional Dependencies: {self.format_fds(self.fds)}\n\n")
        
        # 2NF
        self.output_text.insert(tk.END, "STEP 2: SECOND NORMAL FORM (2NF)\n")
        self.output_text.insert(tk.END, "A relation is in 2NF if it is in 1NF and no non-prime attribute is dependent on a proper subset of any candidate key.\n\n")
        
        # Check for partial dependencies
        partial_deps = self.find_partial_dependencies()
        
        if not partial_deps:
            self.output_text.insert(tk.END, "No partial dependencies found. The schema is already in 2NF.\n\n")
            schema_2nf = schema_1nf
        else:
            self.output_text.insert(tk.END, "Partial dependencies found:\n")
            for subset, attrs in partial_deps:
                self.output_text.insert(tk.END, f"  {', '.join(subset)} -> {', '.join(attrs)}\n")
                
            # Decompose to remove partial dependencies
            schema_2nf = self.decompose_2nf(schema_1nf, partial_deps)
            
            self.output_text.insert(tk.END, "\n2NF Schema:\n")
            for i, rel in enumerate(schema_2nf, 1):
                self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
            
        # 3NF
        self.output_text.insert(tk.END, "\nSTEP 3: THIRD NORMAL FORM (3NF)\n")
        self.output_text.insert(tk.END, "A relation is in 3NF if it is in 2NF and every non-prime attribute is non-transitively dependent on every candidate key.\n\n")
        
        # Find transitive dependencies in each 2NF relation
        transitive_deps = []
        for relation in schema_2nf:
            rel_trans_deps = self.find_transitive_dependencies(relation)
            if rel_trans_deps:
                transitive_deps.extend(rel_trans_deps)
                
        if not transitive_deps:
            self.output_text.insert(tk.END, "No transitive dependencies found. The schema is already in 3NF.\n\n")
            schema_3nf = schema_2nf
        else:
            self.output_text.insert(tk.END, "Transitive dependencies found:\n")
            for a, b, c in transitive_deps:
                self.output_text.insert(tk.END, f"  {a} -> {b} -> {c}\n")
                
            # Decompose to remove transitive dependencies
            schema_3nf = self.decompose_3nf(schema_2nf, transitive_deps)
            
            self.output_text.insert(tk.END, "\n3NF Schema:\n")
            for i, rel in enumerate(schema_3nf, 1):
                self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
        
        # BCNF
        self.output_text.insert(tk.END, "\nSTEP 4: BOYCE-CODD NORMAL FORM (BCNF)\n")
        self.output_text.insert(tk.END, "A relation is in BCNF if for every non-trivial functional dependency X -> Y, X is a superkey.\n\n")
        
        # Check for BCNF violations in each 3NF relation
        bcnf_violations = []
        for relation in schema_3nf:
            rel_violations = self.find_bcnf_violations(relation)
            if rel_violations:
                bcnf_violations.extend(rel_violations)
                
        if not bcnf_violations:
            self.output_text.insert(tk.END, "No BCNF violations found. The schema is already in BCNF.\n\n")
            schema_bcnf = schema_3nf
        else:
            self.output_text.insert(tk.END, "BCNF violations found:\n")
            for lhs, rhs, rel in bcnf_violations:
                self.output_text.insert(tk.END, f"  In relation ({', '.join(rel)}): {', '.join(lhs)} -> {', '.join(rhs)}\n")
                
            # Decompose to eliminate BCNF violations
            schema_bcnf = self.decompose_bcnf(schema_3nf, bcnf_violations)
            
            self.output_text.insert(tk.END, "\nBCNF Schema:\n")
            for i, rel in enumerate(schema_bcnf, 1):
                self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
        
        # Summary
        self.output_text.insert(tk.END, "\n=== NORMALIZATION SUMMARY ===\n\n")
        self.output_text.insert(tk.END, f"Original Schema: R({', '.join(self.attributes)})\n\n")
        
        self.output_text.insert(tk.END, "1NF Schema:\n")
        for i, rel in enumerate(schema_1nf, 1):
            self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
            
        self.output_text.insert(tk.END, "\n2NF Schema:\n")
        for i, rel in enumerate(schema_2nf, 1):
            self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
            
        self.output_text.insert(tk.END, "\n3NF Schema:\n")
        for i, rel in enumerate(schema_3nf, 1):
            self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
            
        self.output_text.insert(tk.END, "\nBCNF Schema:\n")
        for i, rel in enumerate(schema_bcnf, 1):
            self.output_text.insert(tk.END, f"R{i}({', '.join(rel)})\n")
    
    def find_partial_dependencies(self):
        """Find partial dependencies in the schema."""
        partial_deps = []
        
        # If no candidate keys with 2+ attributes, no partial dependencies
        multi_attr_keys = [key for key in self.candidate_keys if len(key) >= 2]
        if not multi_attr_keys:
            return partial_deps
            
        for key in multi_attr_keys:
            # Generate all proper subsets of the key
            proper_subsets = []
            for r in range(1, len(key)):
                proper_subsets.extend(combinations(key, r))
                
            for subset in proper_subsets:
                subset_set = set(subset)
                closure = self.compute_closure(subset_set)
                
                # Check if this subset determines any non-key attributes
                determined_attrs = closure - subset_set
                non_key_attrs = determined_attrs - set().union(*self.candidate_keys)
                
                if non_key_attrs:
                    partial_deps.append((subset, non_key_attrs))
                    
        return partial_deps
    
    def decompose_2nf(self, schema_1nf, partial_deps):
        """Decompose a 1NF schema into 2NF based on partial dependencies."""
        results = []
        
        for relation in schema_1nf:
            # Check if this relation has any partial dependencies
            rel_partial_deps = []
            for subset, attrs in partial_deps:
                if set(subset).issubset(relation) and attrs.issubset(relation):
                    rel_partial_deps.append((subset, attrs))
                    
            if not rel_partial_deps:
                # This relation has no partial dependencies, keep it as is
                results.append(relation)
            else:
                # This relation has partial dependencies, decompose it
                # Start with subsets that have partial dependencies
                processed_attrs = set()
                
                for subset, attrs in rel_partial_deps:
                    # Create a new relation with the subset and its dependent attributes
                    new_rel = set(subset) | attrs
                    results.append(new_rel)
                    processed_attrs.update(new_rel)
                
                # Create a relation with any remaining attributes plus the key attributes
                remaining_attrs = relation - processed_attrs
                if remaining_attrs:
                    # Add a candidate key to the remaining attributes
                    for key in self.candidate_keys:
                        if key.issubset(relation):
                            remaining_attrs.update(key)
                            break
                    
                    results.append(remaining_attrs)
        
        return results
    
    def find_transitive_dependencies(self, relation):
        """Find transitive dependencies in a relation."""
        transitive_deps = []
        rel_attrs = set(relation)
        
        # Get functional dependencies that apply to this relation
        rel_fds = [(lhs, rhs) for lhs, rhs in self.fds 
                  if lhs.issubset(rel_attrs) and rhs.issubset(rel_attrs)]
        
        # Find candidate keys for this relation
        rel_keys = []
        for attr_set in self.powerset(relation):
            if not attr_set:
                continue
                
            closure = self.compute_closure(attr_set) & rel_attrs
            if closure == rel_attrs:
                # Check if it's minimal
                is_minimal = True
                for key in rel_keys:
                    if key.issubset(attr_set):
                        is_minimal = False
                        break
                
                if is_minimal:
                    rel_keys.append(frozenset(attr_set))
        
        # Find prime attributes (those that are part of any candidate key)
        prime_attrs = set()
        for key in rel_keys:
            prime_attrs.update(key)
            
        # Check for transitive dependencies A -> B -> C where C is not part of a key
        for lhs1, rhs1 in rel_fds:
            for lhs2, rhs2 in rel_fds:
                # Check if rhs1 is the same as lhs2 (A -> B -> C pattern)
                if rhs1.issubset(lhs2) and lhs1.isdisjoint(rhs2):
                    # Check if rhs2 contains non-prime attributes
                    non_prime_rhs = rhs2 - prime_attrs
                    if non_prime_rhs:
                        for attr1 in lhs1:
                            for attr2 in rhs1 & lhs2:
                                for attr3 in non_prime_rhs:
                                    transitive_deps.append((attr1, attr2, attr3))
        
        return transitive_deps
    
    def decompose_3nf(self, schema_2nf, transitive_deps):
        """Decompose a 2NF schema into 3NF based on transitive dependencies."""
        results = []
        
        for relation in schema_2nf:
            # Check if this relation has any transitive dependencies
            rel_trans_deps = []
            for a, b, c in transitive_deps:
                if {a, b, c}.issubset(relation):
                    rel_trans_deps.append((a, b, c))
                    
            if not rel_trans_deps:
                # This relation has no transitive dependencies, keep it as is
                results.append(relation)
            else:
                # This relation has transitive dependencies, decompose it
                processed = set()
                
                for a, b, c in rel_trans_deps:
                    # Create new relations to break the transitive dependency
                    # One with A -> B
                    r1 = {a, b}
                    # One with B -> C
                    r2 = {b, c}
                    
                    if r1 not in results and not any(r1.issubset(r) for r in results):
                        results.append(r1)
                    if r2 not in results and not any(r2.issubset(r) for r in results):
                        results.append(r2)
                        
                    processed.update({a, b, c})
                
                # Add a relation with any remaining attributes
                remaining = relation - processed
                if remaining:
                    # Include a key in the remaining relation
                    key_needed = True
                    for key in self.candidate_keys:
                        if key.issubset(relation) and key.issubset(remaining):
                            key_needed = False
                            break
                            
                    if key_needed:
                        # Add a minimal set of attributes to form a key
                        for key in self.candidate_keys:
                            if key.issubset(relation):
                                remaining.update(key)
                                break
                                
                    results.append(remaining)
        
        return results
    
    def find_bcnf_violations(self, relation):
        """Find BCNF violations in a relation."""
        violations = []
        rel_attrs = set(relation)
        
        # Get functional dependencies that apply to this relation
        rel_fds = [(lhs, rhs) for lhs, rhs in self.fds 
                  if lhs.issubset(rel_attrs) and rhs.issubset(rel_attrs)]
        
        # Find all keys for this relation
        rel_keys = []
        for attr_set in self.powerset(relation):
            if not attr_set:
                continue
                
            closure = self.compute_closure(attr_set) & rel_attrs
            if closure == rel_attrs:
                # Check if it's minimal
                is_minimal = True
                for key in rel_keys:
                    if key.issubset(attr_set):
                        is_minimal = False
                        break
                
                if is_minimal:
                    rel_keys.append(frozenset(attr_set))
        
        # Check each FD to see if it violates BCNF
        for lhs, rhs in rel_fds:
            # Skip trivial FDs
            if rhs.issubset(lhs):
                continue
                
            # Check if LHS is a superkey for the relation
            is_superkey = False
            for key in rel_keys:
                if key.issubset(lhs):
                    is_superkey = True
                    break
                    
            if not is_superkey:
                violations.append((lhs, rhs, relation))
                
        return violations
    
    def decompose_bcnf(self, schema_3nf, bcnf_violations):
        """Decompose a 3NF schema into BCNF based on violations."""
        results = list(schema_3nf)  # Start with 3NF schema
        
        # Process each violation
        for lhs, rhs, relation in bcnf_violations:
            # Find which relation in the results contains this violation
            for i, rel in enumerate(results):
                if set(relation).issubset(rel):
                    # Split the relation based on the violation
                    r1 = lhs | rhs  # Relation with the FD
                    r2 = (lhs | (set(rel) - rhs))  # Relation with the key
                    
                    # Replace the original relation
                    results[i] = r1
                    results.append(r2)
                    break
        
        # Remove any redundant relations (completely contained in others)
        i = 0
        while i < len(results):
            redundant = False
            for j, rel in enumerate(results):
                if i != j and results[i].issubset(rel):
                    redundant = True
                    break
            
            if redundant:
                results.pop(i)
            else:
                i += 1
                
        return results
    
    def powerset(self, iterable):
        """Return a generator for all possible subsets of the input iterable."""
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))
    
    def format_fds(self, fds):
        """Format functional dependencies for display."""
        result = []
        for lhs, rhs in fds:
            result.append(f"{', '.join(sorted(lhs))} -> {', '.join(sorted(rhs))}")
        return "; ".join(result)
    
    def clear_results(self):
        """Clear the output text area."""
        self.output_text.delete("1.0", tk.END)

def main():
    root = tk.Tk()
    app = NormalizationApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
