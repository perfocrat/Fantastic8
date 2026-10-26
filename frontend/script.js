const geometryDefinitions = {
    rectangle: {
        name: "Rectangle",
        dimensions: [
            { key: "width", label: "Width", value: 8, unit: "nm" },
            { key: "height", label: "Height", value: 5, unit: "nm" }
        ]
    },
    triangle: {
        name: "Equilateral Triangle",
        dimensions: [{ key: "side", label: "Side length", value: 7, unit: "nm" }]
    },
    circle: {
        name: "Circle",
        dimensions: [{ key: "radius", label: "Radius", value: 4, unit: "nm" }]
    },
    sector: {
        name: "Circular Sector",
        dimensions: [
            { key: "radius", label: "Radius", value: 4, unit: "nm" },
            { key: "angle", label: "Sector angle", value: 90, unit: "°" },
            { key: "startAngle", label: "Starting angle", value: 0, unit: "°" }
        ]
    },
    annulus: {
        name: "Annulus",
        dimensions: [
            { key: "innerRadius", label: "Inner radius", value: 2, unit: "nm" },
            { key: "outerRadius", label: "Outer radius", value: 5, unit: "nm" }
        ]
    },
    ellipse: {
        name: "Ellipse",
        dimensions: [
            { key: "semiMajor", label: "Semi-major axis", value: 5, unit: "nm" },
            { key: "semiMinor", label: "Semi-minor axis", value: 3, unit: "nm" }
        ]
    }
};

const state = {
    geometry: "rectangle",
    quantumNumbers: { n: 1, m: 2 },
    dimensions: Object.fromEntries(
        Object.entries(geometryDefinitions).map(([shape, definition]) => [
            shape,
            Object.fromEntries(definition.dimensions.map(({ key, value }) => [key, value]))
        ])
    )
};

const elements = {
    geometryCards: document.querySelectorAll(".geometry-card"),
    geometryOptions: document.querySelectorAll(".geometry-option"),
    brandHome: document.getElementById("brand-home"),
    landing: document.getElementById("geometry-landing"),
    workspace: document.getElementById("workspace"),
    backButton: document.getElementById("back-to-geometries"),
    workspaceTitle: document.getElementById("workspace-title"),
    geometryName: document.getElementById("active-geometry-name"),
    dimensionUnit: document.getElementById("dimension-unit"),
    dimensionControls: document.getElementById("dimension-controls"),
    form: document.getElementById("parameter-form"),
    calculateButton: document.getElementById("calculate-button"),
    calculationStatus: document.getElementById("calculation-status"),
    backendStatus: document.getElementById("backend-status"),
    badge: document.getElementById("calculation-badge"),
    resultsLabel: document.getElementById("results-label"),
    modelStatus: document.getElementById("model-status"),
    modelStatusDetail: document.getElementById("model-status-detail"),
    energy: document.getElementById("energy-value"),
    energyFoot: document.getElementById("energy-foot"),
    nextEnergy: document.getElementById("next-energy-value"),
    nextEnergyFoot: document.getElementById("next-energy-foot"),
    nSlider: document.getElementById("quantum-n"),
    mSlider: document.getElementById("quantum-m"),
    nValue: document.getElementById("n-value"),
    mValue: document.getElementById("m-value"),
    mMinimumLabel: document.getElementById("m-minimum-label"),
    nMaximumLabel: document.getElementById("n-maximum-label"),
    mMaximumLabel: document.getElementById("m-maximum-label"),
    nLabelContext: document.getElementById("n-label-context"),
    mLabelContext: document.getElementById("m-label-context"),
    nSymbol: document.getElementById("n-symbol"),
    mSymbol: document.getElementById("m-symbol"),
    stateIndexLabel: document.getElementById("state-index-label"),
    constraintNote: document.getElementById("constraint-note"),
    wavefunctionCanvas: document.getElementById("wavefunction-canvas"),
    densityCanvas: document.getElementById("density-canvas"),
    potentialCanvas: document.getElementById("potential-canvas"),
    wavefunctionCaption: document.getElementById("wavefunction-caption"),
    densityCaption: document.getElementById("density-caption"),
    potentialCaption: document.getElementById("potential-caption")
};

function isTriangleSelected() {
    return state.geometry === "triangle";
}

function clearResults() {
    [elements.wavefunctionCanvas, elements.densityCanvas, elements.potentialCanvas].forEach((canvas) => {
        const context = canvas.getContext("2d");
        context.clearRect(0, 0, canvas.width, canvas.height);
    });
    elements.energy.innerHTML = "-- <small>meV</small>";
    elements.nextEnergy.innerHTML = "-- <small>meV</small>";
    elements.energyFoot.textContent = "Not calculated";
    elements.nextEnergyFoot.textContent = "Not available";
    elements.resultsLabel.textContent = "WAITING FOR CALCULATION";
    elements.modelStatus.textContent = "Ready";
    elements.modelStatusDetail.textContent = "Select a shape and calculate.";
    elements.wavefunctionCaption.textContent = "Awaiting calculation";
    elements.densityCaption.textContent = "Awaiting calculation";
    elements.potentialCaption.textContent = "Awaiting calculation";
    elements.calculationStatus.textContent = "Inputs changed. Calculate to update results.";
    elements.badge.innerHTML = "<span></span> NOT CALCULATED";
}

function renderDimensions() {
    const definition = geometryDefinitions[state.geometry];
    elements.dimensionUnit.textContent = definition.dimensions.some(({ unit }) => unit === "°")
        ? "LENGTH / ANGLE"
        : "LENGTH · nm";
    elements.dimensionControls.replaceChildren();

    definition.dimensions.forEach(({ key, label, unit }) => {
        const control = document.createElement("div");
        control.className = "dimension-control";

        const inputId = `dimension-${key}`;
        const fieldLabel = document.createElement("label");
        fieldLabel.htmlFor = inputId;
        fieldLabel.textContent = label;

        const inputWrap = document.createElement("div");
        inputWrap.className = "dimension-input-wrap";

        const input = document.createElement("input");
        input.className = "dimension-input";
        input.id = inputId;
        input.name = key;
        input.type = "number";
        input.min = key === "startAngle" ? "-360" : key === "angle" ? "1" : "0.1";
        if (key === "angle" || key === "startAngle") input.max = "360";
        input.step = unit === "°" ? "1" : "0.1";
        input.required = true;
        input.value = state.dimensions[state.geometry][key];
        input.setAttribute("inputmode", "decimal");
        input.addEventListener("input", () => {
            state.dimensions[state.geometry][key] = input.value === "" ? "" : Number(input.value);
            updateDimensionConstraints();
            clearResults();
        });

        const unitLabel = document.createElement("span");
        unitLabel.className = "dimension-unit";
        unitLabel.textContent = unit;

        inputWrap.append(input, unitLabel);
        control.append(fieldLabel, inputWrap);
        elements.dimensionControls.append(control);
    });
}

function updateDimensionConstraints() {
    const inputs = Object.fromEntries(
        [...elements.dimensionControls.querySelectorAll("input")].map((input) => [input.name, input])
    );
    const dimensions = state.dimensions[state.geometry];
    const pairs = state.geometry === "annulus"
        ? ["innerRadius", "outerRadius"]
        : state.geometry === "ellipse"
            ? ["semiMinor", "semiMajor"]
            : null;
    if (!pairs) return;

    const [smallerKey, largerKey] = pairs;
    const smaller = dimensions[smallerKey];
    const larger = dimensions[largerKey];
    if (inputs[smallerKey] && inputs[largerKey]) {
        inputs[smallerKey].max = String(Math.max(0.1, Number(larger) - 0.1));
        inputs[largerKey].min = String(Number(smaller) + 0.1);
        const validPair = Number.isFinite(smaller) && Number.isFinite(larger) && smaller > 0 && smaller < larger;
        inputs[smallerKey].setCustomValidity(validPair ? "" : `${smallerKey} must be smaller than ${largerKey}.`);
        inputs[largerKey].setCustomValidity(validPair ? "" : `${largerKey} must be larger than ${smallerKey}.`);
    }
}

function quantumLabels() {
    const labels = {
        rectangle: { m: "m", n: "n", mText: "width index", nText: "height index", state: "RECTANGULAR STATE" },
        triangle: { m: "m", n: "n", mText: "upper index", nText: "triangle index", state: "TRIANGLE STATE" },
        circle: { m: "m", n: "n", mText: "angular mode", nText: "radial index", state: "CIRCULAR MODE" },
        sector: { m: "ell", n: "s", mText: "angular index", nText: "radial index", state: "SECTOR MODE" },
        annulus: { m: "l", n: "q", mText: "angular momentum", nText: "radial index", state: "ANNULAR MODE" },
        ellipse: { m: "m", n: "j", mText: "angular mode", nText: "radial index", state: "ELLIPTIC MODE" }
    };
    return labels[state.geometry];
}

function renderQuantumNumbers() {
    const triangleSelected = isTriangleSelected();
    const circleLike = ["circle", "annulus", "ellipse"].includes(state.geometry);
    const minimumM = circleLike ? 0 : 1;
    const maximumM = triangleSelected ? 5 : 10;
    state.quantumNumbers.m = Math.min(maximumM, Math.max(minimumM, state.quantumNumbers.m));
    state.quantumNumbers.n = Math.min(triangleSelected ? state.quantumNumbers.m : 10, state.quantumNumbers.n);
    const labels = quantumLabels();

    elements.mSlider.min = String(minimumM);
    elements.mSlider.max = String(maximumM);
    elements.mSlider.value = String(state.quantumNumbers.m);
    elements.nSlider.min = "1";
    elements.nSlider.max = String(triangleSelected ? state.quantumNumbers.m : 10);
    elements.nSlider.value = String(state.quantumNumbers.n);
    elements.nValue.textContent = String(state.quantumNumbers.n);
    elements.mValue.textContent = String(state.quantumNumbers.m);
    elements.nMaximumLabel.textContent = String(triangleSelected ? state.quantumNumbers.m : 10).padStart(2, "0");
    elements.mMinimumLabel.textContent = String(minimumM).padStart(2, "0");
    elements.mMaximumLabel.textContent = String(maximumM).padStart(2, "0");
    elements.nLabelContext.textContent = labels.nText;
    elements.mLabelContext.textContent = labels.mText;
    elements.nSymbol.textContent = labels.n;
    elements.mSymbol.textContent = labels.m;
    elements.stateIndexLabel.textContent = labels.state;
    elements.constraintNote.hidden = !triangleSelected;
}

function renderGeometry() {
    const definition = geometryDefinitions[state.geometry];
    elements.geometryOptions.forEach((option) => {
        const selected = option.dataset.geometry === state.geometry;
        option.classList.toggle("is-active", selected);
        option.setAttribute("aria-pressed", String(selected));
    });
    elements.workspaceTitle.textContent = definition.name;
    elements.geometryName.textContent = definition.name;
    renderDimensions();
    updateDimensionConstraints();
    renderQuantumNumbers();
}

function selectGeometry(geometry) {
    if (!geometryDefinitions[geometry]) return;
    state.geometry = geometry;
    renderGeometry();
    clearResults();
    showWorkspace();
}

function showWorkspace() {
    elements.landing.hidden = true;
    elements.workspace.hidden = false;
    elements.workspaceTitle.focus?.({ preventScroll: true });
    window.scrollTo({ top: 0, behavior: "smooth" });
}

function showLanding() {
    elements.workspace.hidden = true;
    elements.landing.hidden = false;
    window.scrollTo({ top: 0, behavior: "smooth" });
    elements.landing.querySelector(".geometry-card")?.focus({ preventScroll: true });
}

function onMChange(event) {
    state.quantumNumbers.m = Number(event.currentTarget.value);
    renderQuantumNumbers();
    clearResults();
}

function onNChange(event) {
    state.quantumNumbers.n = Number(event.currentTarget.value);
    renderQuantumNumbers();
    clearResults();
}

function buildRequestPayload() {
    const dimensions = { ...state.dimensions[state.geometry] };
    if (Object.values(dimensions).some((value) => value === "" || !Number.isFinite(value))) {
        throw new Error("Enter a valid number for every geometry parameter.");
    }
    return {
        geometry: state.geometry,
        dimensions,
        quantumNumbers: { ...state.quantumNumbers }
    };
}

function formatEnergy(value) {
    return Number(value).toPrecision(5);
}

function drawField(canvas, values, mask, kind) {
    const rows = values.length;
    const columns = values[0]?.length ?? 0;
    if (!rows || !columns || mask.length !== rows || mask.some((row) => row.length !== columns)) {
        throw new Error("The solver returned an invalid visualization grid.");
    }
    canvas.width = columns;
    canvas.height = rows;
    const context = canvas.getContext("2d", { alpha: true });
    const image = context.createImageData(columns, rows);
    let maxValue = 0;
    values.forEach((row, rowIndex) => row.forEach((value, columnIndex) => {
        if (mask[rowIndex][columnIndex] && Number.isFinite(value)) {
            maxValue = Math.max(maxValue, Math.abs(value));
        }
    }));
    if (kind !== "potential" && (!Number.isFinite(maxValue) || maxValue === 0)) {
        throw new Error("The solver returned an empty field.");
    }

    values.forEach((row, rowIndex) => row.forEach((value, columnIndex) => {
        const canvasRow = rows - rowIndex - 1;
        const offset = (canvasRow * columns + columnIndex) * 4;
        if (!mask[rowIndex][columnIndex] || !Number.isFinite(value)) {
            image.data[offset + 3] = 0;
            return;
        }
        const normalized = Math.min(1, Math.abs(value) / maxValue);
        if (kind === "wave") {
            image.data[offset] = value < 0 ? 232 : 45;
            image.data[offset + 1] = value < 0 ? 153 : 215;
            image.data[offset + 2] = value < 0 ? 83 : 190;
            image.data[offset + 3] = Math.round(45 + normalized * 210);
        } else if (kind === "density") {
            const intensity = Math.sqrt(normalized);
            image.data[offset] = Math.round(14 + 214 * intensity);
            image.data[offset + 1] = Math.round(45 + 190 * intensity);
            image.data[offset + 2] = Math.round(70 + 145 * (1 - intensity));
            image.data[offset + 3] = Math.round(35 + intensity * 220);
        } else {
            image.data[offset] = 42;
            image.data[offset + 1] = 206;
            image.data[offset + 2] = 182;
            image.data[offset + 3] = 230;
        }
    }));
    context.putImageData(image, 0, 0);
}

function renderResult(result) {
    const mask = result.coordinates?.mask;
    drawField(elements.wavefunctionCanvas, result.wavefunction, mask, "wave");
    drawField(elements.densityCanvas, result.probability_density, mask, "density");
    drawField(elements.potentialCanvas, result.potential?.values, mask, "potential");

    elements.energy.innerHTML = `${formatEnergy(result.energy_mev)} <small>meV</small>`;
    const labels = quantumLabels();
    const selectedM = result.inputs.m ?? result.inputs.angular_momentum;
    const selectedN = result.inputs.n ?? result.inputs.radial_index;
    elements.energyFoot.textContent = `${labels.m}=${selectedM}, ${labels.n}=${selectedN}`;
    const nextState = result.states
        .filter((item) => !item.selected)
        .sort((first, second) => first.E - second.E)[0];
    if (nextState) {
        elements.nextEnergy.innerHTML = `${formatEnergy(nextState.energy_mev)} <small>meV</small>`;
        elements.nextEnergyFoot.textContent = `${labels.m}=${nextState.m}, ${labels.n}=${nextState.n}`;
    } else {
        elements.nextEnergy.innerHTML = "-- <small>meV</small>";
        elements.nextEnergyFoot.textContent = "No other listed state";
    }
    elements.resultsLabel.textContent = "SOLVER RESULTS · MEASURED IN MEV";
    elements.modelStatus.textContent = "Calculated";
    elements.modelStatusDetail.textContent = `${result.shape} · ${state.quantumNumbers.m}, ${state.quantumNumbers.n}`;
    elements.calculationStatus.textContent = "Calculation complete.";
    elements.badge.innerHTML = "<span></span> CALCULATED";
    elements.wavefunctionCaption.textContent = `${mask[0].length} × ${mask.length} solver grid`;
    elements.densityCaption.textContent = `${mask[0].length} × ${mask.length} solver grid`;
    elements.potentialCaption.textContent = result.potential.outside_domain === "infinite" ? "Infinite outside boundary" : "Domain returned by solver";
}

async function calculate(event) {
    event.preventDefault();
    if (!elements.form.reportValidity()) return;
    let payload;
    try {
        payload = buildRequestPayload();
    } catch (error) {
        elements.calculationStatus.textContent = error.message;
        return;
    }

    elements.calculateButton.disabled = true;
    elements.calculateButton.innerHTML = "Calculating…";
    elements.calculationStatus.textContent = "Running the selected shape solver…";
    elements.modelStatus.textContent = "Calculating";
    elements.modelStatusDetail.textContent = "Please wait for the numerical solution.";
    elements.badge.innerHTML = "<span></span> CALCULATING";
    try {
        const response = await fetch("/api/calculate", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload)
        });
        const result = await response.json();
        if (!response.ok) throw new Error(result.error || "The calculation request failed.");
        renderResult(result);
    } catch (error) {
        elements.calculationStatus.textContent = error.message || "Unable to reach the physics solver.";
        elements.modelStatus.textContent = "Calculation failed";
        elements.modelStatusDetail.textContent = "Correct the inputs or check the local server.";
        elements.badge.innerHTML = "<span></span> ERROR";
    } finally {
        elements.calculateButton.disabled = false;
        elements.calculateButton.innerHTML = "<span aria-hidden=\"true\">↗</span> Calculate state";
    }
}

async function checkBackend() {
    try {
        const response = await fetch("/api/health");
        if (!response.ok) throw new Error("Backend unavailable");
        elements.backendStatus.innerHTML = "<span></span> Solver connected";
    } catch {
        elements.backendStatus.innerHTML = "<span></span> Solver unavailable";
        elements.calculationStatus.textContent = "Start the local server to calculate states.";
    }
}

elements.geometryOptions.forEach((option) => {
    option.addEventListener("click", () => selectGeometry(option.dataset.geometry));
});
elements.geometryCards.forEach((card) => {
    card.addEventListener("click", () => selectGeometry(card.dataset.geometry));
});
elements.backButton.addEventListener("click", showLanding);
elements.brandHome.addEventListener("click", (event) => {
    if (!elements.workspace.hidden) {
        event.preventDefault();
        showLanding();
    }
});
elements.mSlider.addEventListener("input", onMChange);
elements.nSlider.addEventListener("input", onNChange);
elements.form.addEventListener("submit", calculate);

renderGeometry();
checkBackend();