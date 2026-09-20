(() => {
    const page = document.querySelector('#articlePost.page-teaching');
    if (!page) return;

    function labelCount(count, singular, plural = `${singular}s`) {
        return `${count.toLocaleString('en-US')} ${count === 1 ? singular : plural}`;
    }

    function appendCount(heading, text) {
        const count = document.createElement('span');
        count.className = 'teaching-count';
        count.textContent = `(${text})`;
        heading.append(count);
    }

    page.querySelectorAll('[data-student-list]').forEach(heading => {
        const list = page.querySelector(heading.dataset.studentList);
        if (!list) return;
        const count = list.querySelectorAll('.student-card, .graduate-item').length;
        appendCount(heading.querySelector('span'), labelCount(count, 'student'));
    });

    function summarizeAwards(group, heading) {
        // Count each certificate or scholarship year, including repeat winners.
        const awards = group.querySelectorAll('.cert-item, .scholarship-recipients time').length;
        appendCount(heading, labelCount(awards, 'award'));
    }

    page.querySelectorAll('#supervised-students ~ .award-category').forEach(group => {
        const heading = group.querySelector('.award-title');
        if (heading) summarizeAwards(group, heading);
        group.querySelectorAll('.award-sublist > li').forEach(item => {
            const title = item.querySelector('.scholarship-title');
            if (title) summarizeAwards(item, title);
        });
    });

    const siteNav = document.querySelector('#navBarTop');
    function updateAnchorOffset() {
        const height = siteNav ? siteNav.getBoundingClientRect().height : 0;
        page.style.setProperty('--teaching-anchor-offset', `${height + 20}px`);
    }
    updateAnchorOffset();
    if (siteNav && 'ResizeObserver' in window) new ResizeObserver(updateAnchorOffset).observe(siteNav);
    window.addEventListener('resize', updateAnchorOffset);
})();
