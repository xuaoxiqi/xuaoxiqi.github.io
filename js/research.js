(() => {
    const page = document.querySelector('#articlePost.page-research');
    if (!page) return;

    const conferences = page.querySelector('[data-conf-section="conferences"]');
    const stats = page.querySelector('[data-conf-stats]');
    function updateResearchStats() {
        const overviewCounts = {
            'count-research-projects': page.querySelectorAll('.project-card').length,
            'count-research-conferences': conferences ? conferences.children.length : 0,
            'count-research-seminars': page.querySelectorAll('[data-conf-section="seminars"] > li').length
        };
        Object.entries(overviewCounts).forEach(([id, count]) => {
            const value = document.getElementById(id);
            if (value) value.textContent = count;
        });

        if (!conferences || !stats) return;
        const featured = [...conferences.children].filter(item =>
            item.querySelector('.badge-invite, .badge-keynote')
        );
        stats.querySelector('[data-stat="featured"]').textContent = featured.length;
        stats.hidden = false;
    }
    updateResearchStats();
    if (conferences) {
        new MutationObserver(updateResearchStats).observe(conferences, {
            childList: true,
            subtree: true,
            attributes: true,
            attributeFilter: ['class']
        });
    }

    page.querySelectorAll('[data-start]').forEach(project => {
        if (new Date() < new Date(`${project.dataset.start}T00:00:00`)) return;
        const badge = project.querySelector('.status-badge');
        badge.textContent = 'Active';
        badge.classList.replace('badge-upcoming', 'badge-active');
    });

})();
